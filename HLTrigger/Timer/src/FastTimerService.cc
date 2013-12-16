// FIXME
// we are by-passing the ME's when filling the plots, so we might need to call the ME's update() by hand


// system headers
#ifdef __linux
#include <time.h>
#else
typedef int clockid_t;
#define CLOCK_REALTIME               0
#define CLOCK_MONOTONIC              1
#define CLOCK_PROCESS_CPUTIME_ID     2
#define CLOCK_THREAD_CPUTIME_ID      3
#endif

// C++ headers
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

// boost headers
#include <boost/format.hpp>

// CMSSW headers
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/ModuleCallingContext.h"
#include "FWCore/ServiceRegistry/interface/PathContext.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ServiceRegistry/interface/StreamContext.h"
#include "FWCore/ServiceRegistry/interface/SystemBounds.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "HLTrigger/Timer/interface/FastTimerService.h"


// file-static methods to fill a vector of strings with "(dup.) (...)" entries
static
void fill_dups(std::vector<std::string> & dups, unsigned int size) {
  dups.reserve(size);
  for (unsigned int i = dups.size(); i < size; ++i)
    dups.push_back( (boost::format("(dup.) (%d)") % i).str() );
}


FastTimerService::FastTimerService(const edm::ParameterSet & config, edm::ActivityRegistry & registry) :
  // configuration
  m_use_realtime(                config.getUntrackedParameter<bool>(     "useRealTimeClock"         ) ),
  m_enable_timing_paths(         config.getUntrackedParameter<bool>(     "enableTimingPaths"        ) ),
  m_enable_timing_modules(       config.getUntrackedParameter<bool>(     "enableTimingModules"      ) ),
  m_enable_timing_exclusive(     config.getUntrackedParameter<bool>(     "enableTimingExclusive"    ) ),
  m_enable_timing_summary(       config.getUntrackedParameter<bool>(     "enableTimingSummary"      ) ),
  m_skip_first_path(             config.getUntrackedParameter<bool>(     "skipFirstPath"            ) ),
  // dqm configuration
  m_enable_dqm(                  config.getUntrackedParameter<bool>(     "enableDQM"                ) ),
  m_enable_dqm_bypath_active(    config.getUntrackedParameter<bool>(     "enableDQMbyPathActive"    ) ),
  m_enable_dqm_bypath_total(     config.getUntrackedParameter<bool>(     "enableDQMbyPathTotal"     ) ),
  m_enable_dqm_bypath_overhead(  config.getUntrackedParameter<bool>(     "enableDQMbyPathOverhead"  ) ),
  m_enable_dqm_bypath_details(   config.getUntrackedParameter<bool>(     "enableDQMbyPathDetails"   ) ),
  m_enable_dqm_bypath_counters(  config.getUntrackedParameter<bool>(     "enableDQMbyPathCounters"  ) ),
  m_enable_dqm_bypath_exclusive( config.getUntrackedParameter<bool>(     "enableDQMbyPathExclusive" ) ),
  m_enable_dqm_bymodule(         config.getUntrackedParameter<bool>(     "enableDQMbyModule"        ) ),
  m_enable_dqm_bymoduletype(     config.getUntrackedParameter<bool>(     "enableDQMbyModuleType"    ) ),
  m_enable_dqm_summary(          config.getUntrackedParameter<bool>(     "enableDQMSummary"         ) ),
  m_enable_dqm_byluminosity(     config.getUntrackedParameter<bool>(     "enableDQMbyLuminosity"    ) ),
  m_enable_dqm_byls(             config.getUntrackedParameter<bool>(     "enableDQMbyLumiSection"   ) ),
  m_enable_dqm_bynproc(          config.getUntrackedParameter<bool>(     "enableDQMbyProcesses"     ) ),
  m_nproc_enabled(               false ),
  m_concurrent_runs(             0 ),
  m_concurrent_streams(          0 ),
  m_concurrent_threads(          0 ),
  m_dqm_eventtime_range(         config.getUntrackedParameter<double>(   "dqmTimeRange"             ) ),            // ms
  m_dqm_eventtime_resolution(    config.getUntrackedParameter<double>(   "dqmTimeResolution"        ) ),            // ms
  m_dqm_pathtime_range(          config.getUntrackedParameter<double>(   "dqmPathTimeRange"         ) ),            // ms
  m_dqm_pathtime_resolution(     config.getUntrackedParameter<double>(   "dqmPathTimeResolution"    ) ),            // ms
  m_dqm_moduletime_range(        config.getUntrackedParameter<double>(   "dqmModuleTimeRange"       ) ),            // ms
  m_dqm_moduletime_resolution(   config.getUntrackedParameter<double>(   "dqmModuleTimeResolution"  ) ),            // ms
  m_dqm_luminosity_range(        config.getUntrackedParameter<double>(   "dqmLuminosityRange"       ) / 1.e30),     // cm-2 s-1
  m_dqm_luminosity_resolution(   config.getUntrackedParameter<double>(   "dqmLuminosityResolution"  ) / 1.e30),     // cm-2 s-1
  m_dqm_ls_range(                config.getUntrackedParameter<uint32_t>( "dqmLumiSectionsRange"     ) ),            // lumisections
  m_dqm_path(                    config.getUntrackedParameter<std::string>("dqmPath" ) ),
  m_luminosity_label(            config.getUntrackedParameter<edm::InputTag>("luminosityProduct") ),                // SCAL unpacker
  m_supported_processes(         config.getUntrackedParameter<std::vector<unsigned int> >("supportedProcesses") ),  // 8, 24, 32
  // caching
  m_first_path(),
  m_last_path(),
  m_first_endpath(),
  m_last_endpath(),
  m_is_first_module(false),
  m_is_first_event(true),
  // per-event accounting
  m_event(),
  // per-run and per-job summaries
  m_run_summary(),
  m_job_summary(),
  // DQM - these are initialized at preGlobalBeginRun(), to make sure the DQM service has been loaded
  m_stream()
{
  // enable timers if required by DQM plots
  m_enable_timing_paths     = m_enable_timing_paths         or
                              m_enable_dqm_bypath_active    or
                              m_enable_dqm_bypath_total     or
                              m_enable_dqm_bypath_overhead  or
                              m_enable_dqm_bypath_details   or
                              m_enable_dqm_bypath_counters  or
                              m_enable_dqm_bypath_exclusive;

  m_enable_timing_modules   = m_enable_timing_modules       or
                              m_enable_dqm_bymodule         or
                              m_enable_dqm_bymoduletype     or
                              m_enable_dqm_bypath_total     or
                              m_enable_dqm_bypath_overhead  or
                              m_enable_dqm_bypath_details   or
                              m_enable_dqm_bypath_counters  or
                              m_enable_dqm_bypath_exclusive;

  m_enable_timing_exclusive = m_enable_timing_exclusive     or
                              m_enable_dqm_bypath_exclusive;

  registry.watchPreallocate(       this, & FastTimerService::preallocate );
  registry.watchPreModuleBeginJob( this, & FastTimerService::preModuleBeginJob );
  registry.watchPreGlobalBeginRun( this, & FastTimerService::preGlobalBeginRun );
  registry.watchPostEndJob(        this, & FastTimerService::postEndJob );
  registry.watchPostGlobalEndRun(  this, & FastTimerService::postGlobalEndRun );
  registry.watchPreSourceEvent(    this, & FastTimerService::preSourceEvent );
  registry.watchPostSourceEvent(   this, & FastTimerService::postSourceEvent );
  registry.watchPreEvent(          this, & FastTimerService::preEvent );
  registry.watchPostEvent(         this, & FastTimerService::postEvent );
  // watch per-path events
  registry.watchPrePathEvent(      this, & FastTimerService::prePathEvent );
  registry.watchPostPathEvent(     this, & FastTimerService::postPathEvent );
  // watch per-module events if enabled
  if (m_enable_timing_modules) {
    registry.watchPreModuleEvent(            this, & FastTimerService::preModuleEvent );
    registry.watchPostModuleEvent(           this, & FastTimerService::postModuleEvent );
    registry.watchPreModuleEventDelayedGet(  this, & FastTimerService::preModuleEventDelayedGet );
    registry.watchPostModuleEventDelayedGet( this, & FastTimerService::postModuleEventDelayedGet );
  }
}

FastTimerService::~FastTimerService()
{
}

void FastTimerService::preGlobalBeginRun( edm::GlobalContext const & )
{
  edm::service::TriggerNamesService & tns = * edm::Service<edm::service::TriggerNamesService>();

  // cache the names of the first and last path and endpath
  if (not m_skip_first_path and not tns.getTrigPaths().empty()) {
    m_first_path = tns.getTrigPaths().front();
    m_last_path  = tns.getTrigPaths().back();
  } else if (m_skip_first_path and tns.getTrigPaths().size() > 1) {
    m_first_path = tns.getTrigPaths().at(1);
    m_last_path  = tns.getTrigPaths().back();
  }
  if (not tns.getEndPaths().empty()) {
    m_first_endpath = tns.getEndPaths().front();
    m_last_endpath  = tns.getEndPaths().back();
  }

  uint32_t size_p = tns.getTrigPaths().size();
  uint32_t size_e = tns.getEndPaths().size();
  uint32_t size = size_p + size_e;
  for (uint32_t i = 0; i < size_p; ++i) {
    std::string const & label = tns.getTrigPath(i);
    for (auto & stream: m_stream)
      stream.paths[label].index = i;
  }
  for (uint32_t i = 0; i < size_e; ++i) {
    std::string const & label = tns.getEndPath(i);
    for (auto & stream: m_stream)
      stream.paths[label].index = size_p + i;
  }
  for (auto & timing: m_event)
    timing.paths_interpaths.assign(size + 1, 0);
  for (auto & timing: m_run_summary)
    timing.paths_interpaths.assign(size + 1, 0);
  m_job_summary.paths_interpaths.assign(size + 1, 0);

  // associate to each path all the modules it contains
  for (uint32_t i = 0; i < tns.getTrigPaths().size(); ++i)
    fillPathMap( tns.getTrigPath(i), tns.getTrigPathModules(i) );
  for (uint32_t i = 0; i < tns.getEndPaths().size(); ++i)
    fillPathMap( tns.getEndPath(i), tns.getEndPathModules(i) );


  if (m_enable_dqm) {
    // load the DQM store
    DQMStore * m_dqms = edm::Service<DQMStore>().operator->();

    if (m_dqms == nullptr) {
      // the DQMStore is not available, disable all DQM plots
      m_enable_dqm = false;
      return;
    }

    int eventbins  = (int) std::ceil(m_dqm_eventtime_range  / m_dqm_eventtime_resolution);
    int pathbins   = (int) std::ceil(m_dqm_pathtime_range   / m_dqm_pathtime_resolution);
    int modulebins = (int) std::ceil(m_dqm_moduletime_range / m_dqm_moduletime_resolution);
    int lumibins   = (int) std::ceil(m_dqm_luminosity_range / m_dqm_luminosity_resolution);

    // book MonitorElement's for each stream
    for (unsigned int sid = 0; sid < m_concurrent_streams; ++sid) {
      auto & stream = m_stream[sid];
      // FIXME remove when using threadsafe DQM
      std::string spath = (m_concurrent_streams > 1) ? (boost::format("/stream %02d") % sid).str() : std::string();

      // event summary plots
      if (m_enable_dqm_summary) {
        m_dqms->setCurrentFolder((m_dqm_path + spath));
        stream.dqm.event         = m_dqms->book1D("event",        "Event processing time",         eventbins,  0., m_dqm_eventtime_range)->getTH1F();
        stream.dqm.event         ->StatOverflows(true);
        stream.dqm.event         ->SetXTitle("processing time [ms]");
        stream.dqm.presource     = m_dqms->book1D("presource",    "Pre-Source processing time",    modulebins, 0., m_dqm_moduletime_range)->getTH1F();
        stream.dqm.presource     ->StatOverflows(true);
        stream.dqm.presource     ->SetXTitle("processing time [ms]");
        stream.dqm.source        = m_dqms->book1D("source",       "Source processing time",        modulebins, 0., m_dqm_moduletime_range)->getTH1F();
        stream.dqm.source        ->StatOverflows(true);
        stream.dqm.source        ->SetXTitle("processing time [ms]");
        stream.dqm.postsource    = m_dqms->book1D("postsource",   "Post-Source processing time",   modulebins, 0., m_dqm_moduletime_range)->getTH1F();
        stream.dqm.postsource    ->StatOverflows(true);
        stream.dqm.postsource    ->SetXTitle("processing time [ms]");
        stream.dqm.all_paths     = m_dqms->book1D("all_paths",    "Paths processing time",         eventbins,  0., m_dqm_eventtime_range)->getTH1F();
        stream.dqm.all_paths     ->StatOverflows(true);
        stream.dqm.all_paths     ->SetXTitle("processing time [ms]");
        stream.dqm.all_endpaths  = m_dqms->book1D("all_endpaths", "EndPaths processing time",      pathbins,   0., m_dqm_pathtime_range)->getTH1F();
        stream.dqm.all_endpaths  ->StatOverflows(true);
        stream.dqm.all_endpaths  ->SetXTitle("processing time [ms]");
        stream.dqm.interpaths    = m_dqms->book1D("interpaths",   "Time spent between paths",      pathbins,   0., m_dqm_eventtime_range)->getTH1F();
        stream.dqm.interpaths    ->StatOverflows(true);
        stream.dqm.interpaths    ->SetXTitle("processing time [ms]");
      }

      // event summary plots - summed over nodes with the same number of processes
      if (m_enable_dqm_summary and m_enable_dqm_bynproc) {
        TH1F * plot;
        for (unsigned int n : m_supported_processes) {
          m_dqms->setCurrentFolder( (boost::format("%s/Running %d processes") % (m_dqm_path + spath) % n).str() );
          plot = m_dqms->book1D("event",        "Event processing time",       eventbins,  0., m_dqm_eventtime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.event = plot;
          plot = m_dqms->book1D("presource",    "Pre-Source processing time",  modulebins, 0., m_dqm_moduletime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.presource = plot;
          plot = m_dqms->book1D("source",       "Source processing time",      modulebins, 0., m_dqm_moduletime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.source = plot;
          plot = m_dqms->book1D("postsource",   "Post-Source processing time", modulebins, 0., m_dqm_moduletime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.postsource = plot;
          plot = m_dqms->book1D("all_paths",    "Paths processing time",       eventbins,  0., m_dqm_eventtime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.all_paths = plot;
          plot = m_dqms->book1D("all_endpaths", "EndPaths processing time",    pathbins,   0., m_dqm_pathtime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.all_endpaths = plot;
          plot = m_dqms->book1D("interpaths",   "Time spent between paths",    pathbins,   0., m_dqm_eventtime_range)->getTH1F();
          plot->StatOverflows(true);
          plot->SetXTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc.interpaths = plot;
        }
      }

      // plots by path
      m_dqms->setCurrentFolder((m_dqm_path + spath));
      stream.dqm_paths_active_time     = m_dqms->bookProfile("paths_active_time",    "Additional time spent in each path", size, -0.5, size-0.5, pathbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
      stream.dqm_paths_active_time     ->StatOverflows(true);
      stream.dqm_paths_active_time     ->SetYTitle("processing time [ms]");
      stream.dqm_paths_total_time      = m_dqms->bookProfile("paths_total_time",     "Total time spent in each path",      size, -0.5, size-0.5, pathbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
      stream.dqm_paths_total_time      ->StatOverflows(true);
      stream.dqm_paths_total_time      ->SetYTitle("processing time [ms]");
      stream.dqm_paths_exclusive_time  = m_dqms->bookProfile("paths_exclusive_time", "Exclusive time spent in each path",  size, -0.5, size-0.5, pathbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
      stream.dqm_paths_exclusive_time  ->StatOverflows(true);
      stream.dqm_paths_exclusive_time  ->SetYTitle("processing time [ms]");
      stream.dqm_paths_interpaths      = m_dqms->bookProfile("paths_interpaths",     "Time spent between each path",       size, -0.5, size-0.5, pathbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
      stream.dqm_paths_interpaths      ->StatOverflows(true);
      stream.dqm_paths_interpaths      ->SetYTitle("processing time [ms]");

      for (uint32_t i = 0; i < size_p; ++i) {
        std::string const & label = tns.getTrigPath(i);
        stream.dqm_paths_active_time    ->GetXaxis()->SetBinLabel(i + 1, label.c_str());
        stream.dqm_paths_total_time     ->GetXaxis()->SetBinLabel(i + 1, label.c_str());
        stream.dqm_paths_exclusive_time ->GetXaxis()->SetBinLabel(i + 1, label.c_str());
        stream.dqm_paths_interpaths     ->GetXaxis()->SetBinLabel(i + 1, label.c_str());
      }
      for (uint32_t i = 0; i < size_e; ++i) {
        std::string const & label = tns.getEndPath(i);
        stream.dqm_paths_active_time    ->GetXaxis()->SetBinLabel(i + size_p + 1, label.c_str());
        stream.dqm_paths_total_time     ->GetXaxis()->SetBinLabel(i + size_p + 1, label.c_str());
        stream.dqm_paths_exclusive_time ->GetXaxis()->SetBinLabel(i + size_p + 1, label.c_str());
        stream.dqm_paths_interpaths     ->GetXaxis()->SetBinLabel(i + size_p + 1, label.c_str());
      }

      // per-lumisection plots
      if (m_enable_dqm_byls) {
        m_dqms->setCurrentFolder((m_dqm_path + spath));
        stream.dqm_byls.event        = m_dqms->bookProfile("event_byls",        "Event processing time, by LumiSection",       m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.event        ->StatOverflows(true);
        stream.dqm_byls.event        ->SetXTitle("lumisection");
        stream.dqm_byls.event        ->SetYTitle("processing time [ms]");
        stream.dqm_byls.presource    = m_dqms->bookProfile("presource_byls",    "Pre-Source processing time, by LumiSection",  m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.presource    ->StatOverflows(true);
        stream.dqm_byls.presource    ->SetXTitle("lumisection");
        stream.dqm_byls.presource    ->SetYTitle("processing time [ms]");
        stream.dqm_byls.source       = m_dqms->bookProfile("source_byls",       "Source processing time, by LumiSection",      m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.source       ->StatOverflows(true);
        stream.dqm_byls.source       ->SetXTitle("lumisection");
        stream.dqm_byls.source       ->SetYTitle("processing time [ms]");
        stream.dqm_byls.postsource   = m_dqms->bookProfile("postsource_byls",   "Post-Source processing time, by LumiSection", m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.postsource   ->StatOverflows(true);
        stream.dqm_byls.postsource   ->SetXTitle("lumisection");
        stream.dqm_byls.postsource   ->SetYTitle("processing time [ms]");
        stream.dqm_byls.all_paths    = m_dqms->bookProfile("all_paths_byls",    "Paths processing time, by LumiSection",       m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.all_paths    ->StatOverflows(true);
        stream.dqm_byls.all_paths    ->SetXTitle("lumisection");
        stream.dqm_byls.all_paths    ->SetYTitle("processing time [ms]");
        stream.dqm_byls.all_endpaths = m_dqms->bookProfile("all_endpaths_byls", "EndPaths processing time, by LumiSection",    m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.all_endpaths ->StatOverflows(true);
        stream.dqm_byls.all_endpaths ->SetXTitle("lumisection");
        stream.dqm_byls.all_endpaths ->SetYTitle("processing time [ms]");
        stream.dqm_byls.interpaths   = m_dqms->bookProfile("interpaths_byls",   "Time spent between paths, by LumiSection",    m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byls.interpaths   ->StatOverflows(true);
        stream.dqm_byls.interpaths   ->SetXTitle("lumisection");
        stream.dqm_byls.interpaths   ->SetYTitle("processing time [ms]");
      }

      // per-lumisection plots
      if (m_enable_dqm_byls and m_enable_dqm_bynproc) {
        TProfile * plot;
        for (unsigned int n : m_supported_processes) {
          m_dqms->setCurrentFolder( (boost::format("%s/Running %d processes") % (m_dqm_path + spath) % n).str() );
          plot = m_dqms->bookProfile("event_byls",        "Event processing time, by LumiSection",    m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.event = plot;
          plot = m_dqms->bookProfile("presource_byls",    "Pre-Source processing time, by LumiSection",   m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.presource = plot;
          plot = m_dqms->bookProfile("source_byls",       "Source processing time, by LumiSection",   m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.source = plot;
          plot = m_dqms->bookProfile("postsource_byls",   "Post-Source processing time, by LumiSection",   m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.postsource = plot;
          plot = m_dqms->bookProfile("all_paths_byls",    "Paths processing time, by LumiSection",    m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.all_paths = plot;
          plot = m_dqms->bookProfile("all_endpaths_byls", "EndPaths processing time, by LumiSection", m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.all_endpaths = plot;
          plot = m_dqms->bookProfile("interpaths_byls",   "Time spent between paths, by LumiSection", m_dqm_ls_range, 0.5, m_dqm_ls_range+0.5, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetXTitle("lumisection");
          plot->SetYTitle("processing time [ms]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byls.interpaths = plot;
        }
      }

      // plots vs. instantaneous luminosity
      if (m_enable_dqm_byluminosity) {
        m_dqms->setCurrentFolder((m_dqm_path + spath));
        stream.dqm_byluminosity.event        = m_dqms->bookProfile("event_byluminosity",        "Event processing time vs. instantaneous luminosity",        lumibins, 0., m_dqm_luminosity_range, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.event        ->StatOverflows(true);
        stream.dqm_byluminosity.event        ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.event        ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.presource    = m_dqms->bookProfile("presource_byluminosity",    "Pre-Source processing time vs. instantaneous luminosity",   lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.presource    ->StatOverflows(true);
        stream.dqm_byluminosity.presource    ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.presource    ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.source       = m_dqms->bookProfile("source_byluminosity",       "Source processing time vs. instantaneous luminosity",       lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.source       ->StatOverflows(true);
        stream.dqm_byluminosity.source       ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.source       ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.postsource   = m_dqms->bookProfile("postsource_byluminosity",   "Post-Source processing time vs. instantaneous luminosity",  lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.postsource   ->StatOverflows(true);
        stream.dqm_byluminosity.postsource   ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.postsource   ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.all_paths    = m_dqms->bookProfile("all_paths_byluminosity",    "Paths processing time vs. instantaneous luminosity",        lumibins, 0., m_dqm_luminosity_range, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.all_paths    ->StatOverflows(true);
        stream.dqm_byluminosity.all_paths    ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.all_paths    ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.all_endpaths = m_dqms->bookProfile("all_endpaths_byluminosity", "EndPaths processing time vs. instantaneous luminosity",     lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.all_endpaths ->StatOverflows(true);
        stream.dqm_byluminosity.all_endpaths ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.all_endpaths ->SetYTitle("processing time [ms]");
        stream.dqm_byluminosity.interpaths   = m_dqms->bookProfile("interpaths_byluminosity",   "Time spent between paths vs. instantaneous luminosity",     lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
        stream.dqm_byluminosity.interpaths   ->StatOverflows(true);
        stream.dqm_byluminosity.interpaths   ->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
        stream.dqm_byluminosity.interpaths   ->SetYTitle("processing time [ms]");
      }

      // plots vs. instantaneous luminosity
      if (m_enable_dqm_byluminosity and m_enable_dqm_bynproc) {
        TProfile * plot;
        for (unsigned int n : m_supported_processes) {
          m_dqms->setCurrentFolder( (boost::format("%s/Running %d processes") % (m_dqm_path + spath) % n).str() );
          plot = m_dqms->bookProfile("event_byluminosity",        "Event processing time vs. instantaneous luminosity",       lumibins, 0., m_dqm_luminosity_range, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.event = plot;
          plot = m_dqms->bookProfile("presource_byluminosity",    "Pre-Source processing time vs. instantaneous luminosity",  lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.presource = plot;
          plot = m_dqms->bookProfile("source_byluminosity",       "Source processing time vs. instantaneous luminosity",      lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.source = plot;
          plot = m_dqms->bookProfile("postsource_byluminosity",   "Post-Source processing time vs. instantaneous luminosity", lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.postsource = plot;
          plot = m_dqms->bookProfile("all_paths_byluminosity",    "Paths processing time vs. instantaneous luminosity",       lumibins, 0., m_dqm_luminosity_range, eventbins, 0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.all_paths = plot;
          plot = m_dqms->bookProfile("all_endpaths_byluminosity", "EndPaths processing time vs. instantaneous luminosity",    lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.all_endpaths = plot;
          plot = m_dqms->bookProfile("interpaths_byluminosity",   "Time spent between paths vs. instantaneous luminosity",    lumibins, 0., m_dqm_luminosity_range, pathbins,  0., std::numeric_limits<double>::infinity(), " ")->getTProfile();
          plot->StatOverflows(true);
          plot->SetYTitle("processing time [ms]");
          plot->SetXTitle("instantaneous luminosity [10^{30} cm^{-2}s^{-1}]");
          if (m_concurrent_threads == n) stream.dqm_nproc_byluminosity.interpaths = plot;
        }
      }

      // per-path and per-module accounting
      if (m_enable_timing_paths) {
        m_dqms->setCurrentFolder(((m_dqm_path + spath) + "/Paths"));
        for (auto & keyval: stream.paths) {
          std::string const & pathname = keyval.first;
          PathInfo          & pathinfo = keyval.second;

          if (m_enable_dqm_bypath_active) {
            pathinfo.dqm_active       = m_dqms->book1D(pathname + "_active",       pathname + " active time",            pathbins, 0., m_dqm_pathtime_range)->getTH1F();
            pathinfo.dqm_active       ->StatOverflows(true);
            pathinfo.dqm_active       ->SetXTitle("processing time [ms]");
          }

          if (m_enable_dqm_bypath_total) {
            pathinfo.dqm_total        = m_dqms->book1D(pathname + "_total",        pathname + " total time",             pathbins, 0., m_dqm_pathtime_range)->getTH1F();
            pathinfo.dqm_total        ->StatOverflows(true);
            pathinfo.dqm_total        ->SetXTitle("processing time [ms]");
          }

          if (m_enable_dqm_bypath_overhead) {
            pathinfo.dqm_premodules   = m_dqms->book1D(pathname + "_premodules",   pathname + " pre-modules overhead",   modulebins, 0., m_dqm_moduletime_range)->getTH1F();
            pathinfo.dqm_premodules   ->StatOverflows(true);
            pathinfo.dqm_premodules   ->SetXTitle("processing time [ms]");
            pathinfo.dqm_intermodules = m_dqms->book1D(pathname + "_intermodules", pathname + " inter-modules overhead", modulebins, 0., m_dqm_moduletime_range)->getTH1F();
            pathinfo.dqm_intermodules ->StatOverflows(true);
            pathinfo.dqm_intermodules ->SetXTitle("processing time [ms]");
            pathinfo.dqm_postmodules  = m_dqms->book1D(pathname + "_postmodules",  pathname + " post-modules overhead",  modulebins, 0., m_dqm_moduletime_range)->getTH1F();
            pathinfo.dqm_postmodules  ->StatOverflows(true);
            pathinfo.dqm_postmodules  ->SetXTitle("processing time [ms]");
            pathinfo.dqm_overhead     = m_dqms->book1D(pathname + "_overhead",     pathname + " overhead time",          modulebins, 0., m_dqm_moduletime_range)->getTH1F();
            pathinfo.dqm_overhead     ->StatOverflows(true);
            pathinfo.dqm_overhead     ->SetXTitle("processing time [ms]");
          }

          if (m_enable_dqm_bypath_details or m_enable_dqm_bypath_counters) {
            // book histograms for modules-in-paths statistics

            // find histograms X-axis labels
            uint32_t id;
            std::vector<std::string> const & modules = ((id = tns.findTrigPath(pathname)) != tns.getTrigPaths().size()) ? tns.getTrigPathModules(id) :
                                                       ((id = tns.findEndPath(pathname))  != tns.getEndPaths().size())  ? tns.getEndPathModules(id)  :
                                                       std::vector<std::string>();

            static std::vector<std::string> dup;
            if (modules.size() > dup.size())
              fill_dups(dup, modules.size());

            std::vector<const char *> labels(modules.size(), nullptr);
            for (uint32_t i = 0; i < modules.size(); ++i)
              labels[i] = (pathinfo.modules[i]) ? modules[i].c_str() : dup[i].c_str();

            // book counter histograms
            if (m_enable_dqm_bypath_counters) {
              pathinfo.dqm_module_counter = m_dqms->book1D(pathname + "_module_counter", pathname + " module counter", modules.size(), -0.5, modules.size() - 0.5)->getTH1F();
              // find module labels
              for (uint32_t i = 0; i < modules.size(); ++i) {
                pathinfo.dqm_module_counter->GetXaxis()->SetBinLabel( i+1, labels[i] );
              }
            }
            // book detailed timing histograms
            if (m_enable_dqm_bypath_details) {
              pathinfo.dqm_module_active  = m_dqms->book1D(pathname + "_module_active",  pathname + " module active",  modules.size(), -0.5, modules.size() - 0.5)->getTH1F();
              pathinfo.dqm_module_active  ->SetYTitle("cumulative processing time [ms]");
              pathinfo.dqm_module_total   = m_dqms->book1D(pathname + "_module_total",   pathname + " module total",   modules.size(), -0.5, modules.size() - 0.5)->getTH1F();
              pathinfo.dqm_module_total   ->SetYTitle("cumulative processing time [ms]");
              // find module labels
              for (uint32_t i = 0; i < modules.size(); ++i) {
                pathinfo.dqm_module_active ->GetXaxis()->SetBinLabel( i+1, labels[i] );
                pathinfo.dqm_module_total  ->GetXaxis()->SetBinLabel( i+1, labels[i] );
              }
            }
          }

          // book exclusive path time histograms
          if (m_enable_dqm_bypath_exclusive) {
            pathinfo.dqm_exclusive = m_dqms->book1D(pathname + "_exclusive", pathname + " exclusive time", pathbins, 0., m_dqm_pathtime_range)->getTH1F();
            pathinfo.dqm_exclusive ->StatOverflows(true);
            pathinfo.dqm_exclusive ->SetXTitle("processing time [ms]");
          }

        }
      }

      if (m_enable_dqm_bymodule) {
        m_dqms->setCurrentFolder(((m_dqm_path + spath) + "/Modules"));
        for (auto & keyval: stream.modules) {
          std::string const & label  = keyval.first;
          ModuleInfo        & module = keyval.second;
          module.dqm_active = m_dqms->book1D(label, label, modulebins, 0., m_dqm_moduletime_range)->getTH1F();
          module.dqm_active->StatOverflows(true);
          module.dqm_active->SetXTitle("processing time [ms]");
        }
      }

      if (m_enable_dqm_bymoduletype) {
        m_dqms->setCurrentFolder(((m_dqm_path + spath) + "/ModuleTypes"));
        for (auto & keyval: stream.moduletypes) {
          std::string const & label  = keyval.first;
          ModuleInfo        & module = keyval.second;
          module.dqm_active = m_dqms->book1D(label, label, modulebins, 0., m_dqm_moduletime_range)->getTH1F();
          module.dqm_active->StatOverflows(true);
          module.dqm_active->SetXTitle("processing time [ms]");
        }
      }

    }
  }
}


void
FastTimerService::preallocate(edm::service::SystemBounds const& bounds)
{
  m_concurrent_runs    = bounds.maxNumberOfConcurrentRuns();
  m_concurrent_streams = bounds.maxNumberOfStreams();
  m_concurrent_threads = bounds.maxNumberOfThreads();

  m_event.resize(m_concurrent_streams);
  m_run_summary.resize(m_concurrent_runs);
  m_stream.resize(m_concurrent_streams);

  if (m_enable_dqm_bynproc and std::find(m_supported_processes.begin(), m_supported_processes.end(), m_concurrent_threads) != m_supported_processes.end())
    m_nproc_enabled = true;
}


void
FastTimerService::postEndJob()
{
  /*
  if (m_enable_timing_summary) {
    // print a timing sumary for the run
    edm::service::TriggerNamesService & tns = * edm::Service<edm::service::TriggerNamesService>();

    std::ostringstream out;
    out << std::fixed << std::setprecision(6);
    out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ") << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.presource    / (double) m_job_summary.count << "  Pre-Source"    << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.source       / (double) m_job_summary.count << "  Source"        << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.postsource   / (double) m_job_summary.count << "  Post-Source"   << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.event        / (double) m_job_summary.count << "  Event"         << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.all_paths    / (double) m_job_summary.count << "  all Paths"     << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.all_endpaths / (double) m_job_summary.count << "  all EndPaths"  << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_job_summary.interpaths   / (double) m_job_summary.count << "  between paths" << '\n';
    if (m_enable_timing_modules) {
      double modules_total = 0.;
      for (auto & keyval: m_stream.modules)
        modules_total += keyval.second.summary_active;
      out << "FastReport              " << std::right << std::setw(10) << modules_total / (double) m_job_summary.count << "  all Modules"   << '\n';
    }
    out << '\n';
    if (m_enable_timing_paths and not m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active Path" << '\n';
      for (auto const & name: tns.getTrigPaths())
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active / (double) m_job_summary.count << "  "
            << name << '\n';
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active EndPath" << '\n';
      for (auto const & name: tns.getEndPaths())
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active / (double) m_job_summary.count << "  "
            << name << '\n';
    } else if (m_enable_timing_paths and m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active       Pre-     Inter-  Post-mods   Overhead      Total  Path" << '\n';
      for (auto const & name: tns.getTrigPaths()) {
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active        / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_premodules    / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_intermodules  / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_postmodules   / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_overhead      / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_total         / (double) m_job_summary.count << "  "
            << name << '\n';
      }
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active       Pre-     Inter-  Post-mods   Overhead      Total  EndPath" << '\n';
      for (auto const & name: tns.getEndPaths()) {
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active        / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_premodules    / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_intermodules  / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_postmodules   / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_overhead      / (double) m_job_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_total         / (double) m_job_summary.count << "  "
            << name << '\n';
      }
    }
    out << '\n';
    if (m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      for (auto & keyval: m_stream.modules) {
        std::string const & label  = keyval.first;
        ModuleInfo  const & module = keyval.second;
        out << "FastReport              " << std::right << std::setw(10) << module.summary_active  / (double) m_job_summary.count << "  " << label << '\n';
      }
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      for (auto & keyval: m_stream.moduletypes) {
        std::string const & label  = keyval.first;
        ModuleInfo  const & module = keyval.second;
        out << "FastReport              " << std::right << std::setw(10) << module.summary_active  / (double) m_job_summary.count << "  " << label << '\n';
      }
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
    }
    out << '\n';
    edm::LogVerbatim("FastReport") << out.str();
  }
  */
}


void
FastTimerService::postGlobalEndRun(edm::GlobalContext const &)
{
  /*
  if (m_enable_timing_summary) {
    // print a timing sumary for the run
    edm::service::TriggerNamesService & tns = * edm::Service<edm::service::TriggerNamesService>();

    std::ostringstream out;
    out << std::fixed << std::setprecision(6);
    out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ") << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.presource    / (double) m_run_summary.count << "  Pre-Source"    << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.source       / (double) m_run_summary.count << "  Source"        << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.postsource   / (double) m_run_summary.count << "  Post-Source"   << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.event        / (double) m_run_summary.count << "  Event"         << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.all_paths    / (double) m_run_summary.count << "  all Paths"     << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.all_endpaths / (double) m_run_summary.count << "  all EndPaths"  << '\n';
    out << "FastReport              " << std::right << std::setw(10) << m_run_summary.interpaths   / (double) m_run_summary.count << "  between paths" << '\n';
    if (m_enable_timing_modules) {
      double modules_total = 0.;
      for (auto & keyval: m_stream.modules)
        modules_total += keyval.second.summary_active;
      out << "FastReport              " << std::right << std::setw(10) << modules_total / (double) m_run_summary.count << "  all Modules"   << '\n';
    }
    out << '\n';
    if (m_enable_timing_paths and not m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active Path" << '\n';
      for (auto const & name: tns.getTrigPaths())
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active / (double) m_run_summary.count << "  "
            << name << '\n';
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active EndPath" << '\n';
      for (auto const & name: tns.getEndPaths())
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active / (double) m_run_summary.count << "  "
            << name << '\n';
    } else if (m_enable_timing_paths and m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active       Pre-     Inter-  Post-mods   Overhead      Total  Path" << '\n';
      for (auto const & name: tns.getTrigPaths()) {
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active        / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_premodules    / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_intermodules  / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_postmodules   / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_overhead      / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_total         / (double) m_run_summary.count << "  "
            << name << '\n';
      }
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active       Pre-     Inter-  Post-mods   Overhead      Total  EndPath" << '\n';
      for (auto const & name: tns.getEndPaths()) {
        out << "FastReport              "
            << std::right << std::setw(10) << m_stream.paths[name].summary_active        / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_premodules    / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_intermodules  / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_postmodules   / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_overhead      / (double) m_run_summary.count << " "
            << std::right << std::setw(10) << m_stream.paths[name].summary_total         / (double) m_run_summary.count << "  "
            << name << '\n';
      }
    }
    out << '\n';
    if (m_enable_timing_modules) {
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      for (auto & keyval: m_stream.modules) {
        std::string const & label  = keyval.first;
        ModuleInfo  const & module = keyval.second;
        out << "FastReport              " << std::right << std::setw(10) << module.summary_active  / (double) m_run_summary.count << "  " << label << '\n';
      }
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      out << '\n';
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
      for (auto & keyval: m_stream.moduletypes) {
        std::string const & label  = keyval.first;
        ModuleInfo  const & module = keyval.second;
        out << "FastReport              " << std::right << std::setw(10) << module.summary_active  / (double) m_run_summary.count << "  " << label << '\n';
      }
      out << "FastReport " << (m_use_realtime ? "(real time) " : "(CPU time)  ")    << "     Active  Module" << '\n';
    }
    out << '\n';
    edm::LogVerbatim("FastReport") << out.str();
  }
  */
}

void FastTimerService::preModuleBeginJob(edm::ModuleDescription const & module) {
  // allocate a counter for each module and module type
  for (auto & stream: m_stream) {
    stream.fast_modules[& module]     = & stream.modules[module.moduleLabel()];;
    stream.fast_moduletypes[& module] = & stream.moduletypes[module.moduleName()];
  }
}

void FastTimerService::preEvent(edm::StreamContext const & sc) {
  unsigned int sid = sc.streamID().value();

  // new event, reset the per-event counter
  start(m_timer_event);

  // account the time spent after the source
  m_event[sid].postsource = delta(m_timer_source.getStopTime(), m_timer_event.getStartTime());

  // clear the event counters
  m_event[sid].event        = 0;
  m_event[sid].all_paths    = 0;
  m_event[sid].all_endpaths = 0;
  m_event[sid].interpaths   = 0;
  m_event[sid].paths_interpaths.assign(m_stream[sid].paths.size() + 1, 0);
  for (auto & keyval : m_stream[sid].paths) {
    keyval.second.time_active       = 0.;
    keyval.second.time_exclusive    = 0.;
    keyval.second.time_premodules   = 0.;
    keyval.second.time_intermodules = 0.;
    keyval.second.time_postmodules  = 0.;
    keyval.second.time_total        = 0.;
  }
  for (auto & keyval : m_stream[sid].modules) {
    keyval.second.time_active       = 0.;
    keyval.second.run_in_path       = nullptr;
    keyval.second.counter           = 0;
  }
  for (auto & keyval : m_stream[sid].moduletypes) {
    keyval.second.time_active       = 0.;
    keyval.second.run_in_path       = nullptr;
    keyval.second.counter           = 0;
  }

  // copy the start event timestamp as the end of the previous path
  // used by the inter-path overhead measurement
  m_timer_path.setStopTime( m_timer_event.getStartTime() );
}

void FastTimerService::postEvent(edm::StreamContext const & sc) {
  unsigned int sid = sc.streamID();
  unsigned int rid = sc.runIndex();

  // stop the per-event timer, and account event time
  stop(m_timer_event);
  m_event[sid].event = delta(m_timer_event);

  // the last part of inter-path overhead is the time between the end of the last (end)path and the end of the event processing
  double interpaths = delta(m_timer_path.getStopTime(), m_timer_event.getStopTime());
  m_event[sid].interpaths += interpaths;
  m_event[sid].paths_interpaths[m_stream[sid].paths.size()] = interpaths;

  // keep track of the total number of events and add this event's time to the per-run and per-job summary
  m_event[sid].count = 1;
  m_run_summary[rid] += m_event[sid];
  m_job_summary += m_event[sid];

  // elaborate "exclusive" modules
  if (m_enable_timing_exclusive) {
    for (auto & keyval: m_stream[sid].paths) {
      PathInfo & pathinfo = keyval.second;
      pathinfo.time_exclusive = pathinfo.time_overhead;

      for (uint32_t i = 0; i <= pathinfo.last_run; ++i) {
        ModuleInfo * module = pathinfo.modules[i];
        if (module == 0)
          // this is a module occurring more than once in the same path, skip it after the first occurrence
          continue;
        if ((module->run_in_path == & pathinfo) and (module->counter == 1))
          pathinfo.time_exclusive += module->time_active;
      }
    }
  }

  // done processing the first event
  m_is_first_event = false;

  // fill the DQM plots from the internal buffers
  if (not m_enable_dqm)
    return;

  // fill plots for per-event time by path
  if (m_enable_timing_paths) {

    for (auto & keyval: m_stream[sid].paths) {
      PathInfo & pathinfo = keyval.second;

      m_stream[sid].dqm_paths_active_time->Fill(pathinfo.index, pathinfo.time_active * 1000.);
      if (m_enable_dqm_bypath_active)
        pathinfo.dqm_active      ->Fill(pathinfo.time_active          * 1000.);

      m_stream[sid].dqm_paths_exclusive_time->Fill(pathinfo.index, pathinfo.time_exclusive * 1000.);
      if (m_enable_dqm_bypath_exclusive)
        pathinfo.dqm_exclusive   ->Fill(pathinfo.time_exclusive       * 1000.);

      m_stream[sid].dqm_paths_total_time->Fill(pathinfo.index, pathinfo.time_total * 1000.);
      if (m_enable_dqm_bypath_total)
        pathinfo.dqm_total       ->Fill(pathinfo.time_total           * 1000.);

      // fill path overhead histograms
      if (m_enable_dqm_bypath_overhead) {
        pathinfo.dqm_premodules  ->Fill(pathinfo.time_premodules      * 1000.);
        pathinfo.dqm_intermodules->Fill(pathinfo.time_intermodules    * 1000.);
        pathinfo.dqm_postmodules ->Fill(pathinfo.time_postmodules     * 1000.);
        pathinfo.dqm_overhead    ->Fill(pathinfo.time_overhead        * 1000.);
      }

      // fill detailed timing histograms
      if (m_enable_dqm_bypath_details) {
        for (uint32_t i = 0; i <= pathinfo.last_run; ++i) {
          ModuleInfo * module = pathinfo.modules[i];
          // skip duplicate modules
          if (module == nullptr)
            continue;
          // fill the total time for all non-duplicate modules
          pathinfo.dqm_module_total->Fill(i, module->time_active * 1000.);
          // fill the active time only for module that have actually run in this path
          if (module->run_in_path == & pathinfo)
            pathinfo.dqm_module_active->Fill(i, module->time_active * 1000.);
        }
      }

      // fill path counter histograms
      //   - also for duplicate modules, to properly extract rejection information
      //   - fill the N+1th bin for paths accepting the event, so the FastTimerServiceClient can properly measure the last filter efficiency
      if (m_enable_dqm_bypath_counters) {
        for (uint32_t i = 0; i <= pathinfo.last_run; ++i)
          pathinfo.dqm_module_counter->Fill(i);
        if (pathinfo.accept)
          pathinfo.dqm_module_counter->Fill(pathinfo.modules.size());
      }

    }
  }

  // fill plots for per-event time by module
  if (m_enable_dqm_bymodule) {
    for (auto & keyval : m_stream[sid].modules) {
      ModuleInfo & module = keyval.second;
      module.dqm_active->Fill(module.time_active * 1000.);
    }
  }

  // fill plots for per-event time by module type
  if (m_enable_dqm_bymoduletype) {
    for (auto & keyval : m_stream[sid].moduletypes) {
      ModuleInfo & module = keyval.second;
      module.dqm_active->Fill(module.time_active * 1000.);
    }
  }

  // fill the interpath overhead plot
  for (unsigned int i = 0; i <= m_stream[sid].paths.size(); ++i)
    m_stream[sid].dqm_paths_interpaths->Fill(i, m_event[sid].paths_interpaths[i] * 1000.);

  if (m_enable_dqm_summary) {
    m_stream[sid].dqm.fill(m_event[sid]);

    if (m_nproc_enabled)
      m_stream[sid].dqm_nproc.fill(m_event[sid]);
  }

  if (m_enable_dqm_byls) {
    unsigned int ls = sc.eventID().luminosityBlock();
    m_stream[sid].dqm_byls.fill(ls, m_event[sid] );

    if (m_nproc_enabled)
      m_stream[sid].dqm_nproc_byls.fill(ls, m_event[sid]);
  }

  if (m_enable_dqm_byluminosity) {
    float luminosity = 0.;
    /*
     * FIXME - find an alternative approach, the new signals do not pass the actual Event here
    edm::Handle<LumiScalersCollection> h_luminosity;
    if (event.getByLabel(m_luminosity_label, h_luminosity) and not h_luminosity->empty())
      luminosity = h_luminosity->front().instantLumi();   // in units of 1e30 cm-2s-1
    */
    m_stream[sid].dqm_byluminosity.fill(luminosity, m_event[sid]);

    if (m_nproc_enabled)
      m_stream[sid].dqm_nproc_byluminosity.fill(luminosity, m_event[sid]);
  }

}

void FastTimerService::preSourceEvent(edm::StreamID sid) {
  start(m_timer_source);

  // account the time spent before the source
  m_event[sid].presource = (m_is_first_event) ? 0. : delta(m_timer_event.getStopTime(), m_timer_source.getStartTime());

  // clear the event counters
  m_event[sid].source = 0.;
  m_event[sid].postsource = 0.;
}

void FastTimerService::postSourceEvent(edm::StreamID sid) {
  stop(m_timer_source);

  m_event[sid].source = delta(m_timer_source);
}


void FastTimerService::prePathEvent(edm::StreamContext const & sc, edm::PathContext const & pc) {
  std::string const & path = pc.pathName();
  unsigned int sid = sc.streamID();
  auto & stream = m_stream[sid];

  // prepare to measure the time spent between the beginning of the path and the execution of the first module
  m_is_first_module = true;

  PathMap<PathInfo>::iterator keyval = stream.paths.find(path);
  if (keyval != stream.paths.end()) {
    stream.current_path = & keyval->second;
  } else {
    // should never get here
    stream.current_path = 0;
    edm::LogError("FastTimerService") << "FastTimerService::prePathEvent: unexpected path " << path;
  }

  // time each (end)path
  start(m_timer_path);

  if (path == m_first_path) {
    // this is the first path, start the "all paths" counter
    m_timer_paths.setStartTime( m_timer_path.getStartTime() );
  } else if (path == m_first_endpath) {
    // this is the first endpath, start the "all paths" counter
    m_timer_endpaths.setStartTime(m_timer_path.getStartTime());
  }

  // measure the inter-path overhead as the time elapsed since the end of preiovus path
  // (or the beginning of the event, if this is the first path - see preEvent)
  double interpaths = delta(m_timer_path.getStopTime(), m_timer_path.getStartTime());
  m_event[sid].interpaths += interpaths;
  m_event[sid].paths_interpaths[stream.current_path->index] = interpaths;
}


void FastTimerService::postPathEvent(edm::StreamContext const & sc, edm::PathContext const & pc, edm::HLTPathStatus const & status) {
  std::string const & path = pc.pathName();
  unsigned int        sid = sc.streamID().value();
  auto &              stream = m_stream[sid];

  // time each (end)path
  stop(m_timer_path);

  double active = delta(m_timer_path);

  // if enabled, account each (end)path
  if (m_enable_timing_paths) {

    PathInfo & pathinfo = * stream.current_path;
    pathinfo.time_active = active;
    pathinfo.summary_active += active;

    // measure the time spent between the execution of the last module and the end of the path
    if (m_enable_timing_modules) {
      double pre      = 0.;                 // time spent before the first active module
      double inter    = 0.;                 // time spent between active modules
      double post     = 0.;                 // time spent after the last active module
      double overhead = 0.;                 // time spent before, between, or after modules
      double current  = 0.;                 // time spent in modules active in the current path
      double total    = active;             // total per-path time, including modules already run as part of other paths

      // implementation note:
      // "active"   has already measured all the time spent in this path
      // "current"  will be the sum of the time spent inside each module while running this path, so that
      // "overhead" will be active - current
      // "total"    will be active + the sum of the time spent in non-active modules

      uint32_t last_run = status.index();     // index of the last module run in this path
      for (uint32_t i = 0; i <= last_run; ++i) {
        ModuleInfo * module = pathinfo.modules[i];

        if (module == 0)
          // this is a module occurring more than once in the same path, skip it after the first occurrence
          continue;

        ++module->counter;
        if (module->run_in_path == & pathinfo) {
          current += module->time_active;
        } else {
          total   += module->time_active;
        }

      }

      if (m_is_first_module) {
        // no modules were active during this path, account all the time as overhead
        pre      = 0.;
        inter    = 0.;
        post     = active;
        overhead = active;
      } else {
        // extract overhead information
        pre      = delta(m_timer_path.getStartTime(),  m_timer_first_module);
        post     = delta(m_timer_module.getStopTime(), m_timer_path.getStopTime());
        inter    = active - pre - current - post;
        // take care of numeric precision and rounding errors - the timer is less precise than nanosecond resolution
        if (std::abs(inter) < 1e-9)
          inter = 0.;
        overhead = active - current;
        // take care of numeric precision and rounding errors - the timer is less precise than nanosecond resolution
        if (std::abs(overhead) < 1e-9)
          overhead = 0.;
      }

      pathinfo.time_premodules       = pre;
      pathinfo.time_intermodules     = inter;
      pathinfo.time_postmodules      = post;
      pathinfo.time_overhead         = overhead;
      pathinfo.time_total            = total;
      pathinfo.summary_premodules   += pre;
      pathinfo.summary_intermodules += inter;
      pathinfo.summary_postmodules  += post;
      pathinfo.summary_overhead     += overhead;
      pathinfo.summary_total        += total;
      pathinfo.last_run              = status.index();
      pathinfo.accept                = status.accept();
    }
  }

  if (path == m_last_path) {
    // this is the last path, stop and account the "all paths" counter
    m_timer_paths.setStopTime(m_timer_path.getStopTime());
    m_event[sid].all_paths = delta(m_timer_paths);
  } else if (path == m_last_endpath) {
    // this is the last endpath, stop and account the "all endpaths" counter
    m_timer_endpaths.setStopTime(m_timer_path.getStopTime());
    m_event[sid].all_endpaths = delta(m_timer_endpaths);
  }

}

void FastTimerService::preModuleEvent(edm::StreamContext const &, edm::ModuleCallingContext const &) {
  // this is ever called only if m_enable_timing_modules = true
  assert(m_enable_timing_modules);

  // time each module
  start(m_timer_module);

  if (m_is_first_module) {
    m_is_first_module = false;

    // measure the time spent between the beginning of the path and the execution of the first module
    m_timer_first_module = m_timer_module.getStartTime();
  }
}

void FastTimerService::postModuleEvent(edm::StreamContext const & sc, edm::ModuleCallingContext const & mcc) {
  edm::ModuleDescription const * md = mcc.moduleDescription();
  unsigned int sid = sc.streamID().value();
  auto & stream = m_stream[sid];

  // this is ever called only if m_enable_timing_modules = true
  assert(m_enable_timing_modules);

  // time and account each module
  stop(m_timer_module);

  { // if-scope
    ModuleMap<ModuleInfo*>::iterator keyval = stream.fast_modules.find(md);
    if (keyval != stream.fast_modules.end()) {
      double time = delta(m_timer_module);
      ModuleInfo & module = * keyval->second;
      module.run_in_path     = stream.current_path;
      module.time_active     = time;
      module.summary_active += time;
      // plots are filled post event processing
    } else {
      // should never get here
      if (md == nullptr)
        edm::LogError("FastTimerService") << "FastTimerService::postModuleEvent: invalid module";
      else
        edm::LogError("FastTimerService") << "FastTimerService::postModuleEvent: unexpected module " << md->moduleLabel();
    }
  }

  { // if-scope
    ModuleMap<ModuleInfo*>::iterator keyval = stream.fast_moduletypes.find(md);
    if (keyval != stream.fast_moduletypes.end()) {
      double time = delta(m_timer_module);
      ModuleInfo & module = * keyval->second;
      // module.run_in_path is not meaningful here
      module.time_active    += time;
      module.summary_active += time;
      // plots are filled post event processing
    } else {
      // should never get here
      if (md == nullptr)
        edm::LogError("FastTimerService") << "FastTimerService::postModuleEvent: invalid module";
      else
        edm::LogError("FastTimerService") << "FastTimerService::postModuleEvent: unexpected module " << md->moduleLabel();
    }
  }

}

void FastTimerService::preModuleEventDelayedGet(edm::StreamContext const &, edm::ModuleCallingContext const &) {
  // this is ever called only if m_enable_timing_modules = true
  assert(m_enable_timing_modules);

  // XXX
  // check if the signal has interrupted a module or not
  //  - if it is, pause the time for that module, and prepare timing a new module
  //  - otherwise, ignore it
}

void FastTimerService::postModuleEventDelayedGet(edm::StreamContext const &, edm::ModuleCallingContext const &) {
  // this is ever called only if m_enable_timing_modules = true
  assert(m_enable_timing_modules);

  // XXX
  // check if the signal has interrupted a module or not
  //  - if it is, reasume the timing for the original module
  //  - otherwise, ignore it
}

// associate to a path all the modules it contains
void FastTimerService::fillPathMap(std::string const & name, std::vector<std::string> const & modules) {
  for (auto & stream: m_stream) {

    std::vector<ModuleInfo *> & pathmap = stream.paths[name].modules;
    pathmap.clear();
    pathmap.reserve( modules.size() );
    std::unordered_set<ModuleInfo const *> pool;        // keep track of inserted modules
    for (auto const & module: modules) {
      // fix the name of negated or ignored modules
      std::string const & label = (module[0] == '!' or module[0] == '-') ? module.substr(1) : module;

      auto const & it = stream.modules.find(label);
      if (it == stream.modules.end()) {
        // no matching module was found
        pathmap.push_back( 0 );
      } else if (pool.insert(& it->second).second) {
        // new module
        pathmap.push_back(& it->second);
      } else {
        // duplicate module
        pathmap.push_back( 0 );
      }
    }

  }
}


// query the current module/path/event
// Note: these functions incur in a "per-call timer overhead" (see above), currently of the order of 340ns

// return the time spent since the last preModuleEvent() event
double FastTimerService::currentModuleTime(edm::StreamID) const {
  return std::chrono::duration_cast<std::chrono::duration<double>>(m_timer_module.untilNow()).count();
}

// return the time spent since the last prePathEvent() event
double FastTimerService::currentPathTime(edm::StreamID) const {
  return std::chrono::duration_cast<std::chrono::duration<double>>(m_timer_path.untilNow()).count();
}

// return the time spent since the last preEvent() event
double FastTimerService::currentEventTime(edm::StreamID) const {
  return std::chrono::duration_cast<std::chrono::duration<double>>(m_timer_event.untilNow()).count();
}

// query the time spent in a module (available after the module has run)
double FastTimerService::queryModuleTime(edm::StreamID sid, const edm::ModuleDescription & module) const {
  ModuleMap<ModuleInfo *>::const_iterator keyval = m_stream[sid].fast_modules.find(& module);
  if (keyval != m_stream[sid].fast_modules.end()) {
    return keyval->second->time_active;
  } else {
    edm::LogError("FastTimerService") << "FastTimerService::queryModuleTime: unexpected module " << module.moduleLabel();
    return 0.;
  }
}

// query the time spent in a module (available after the module has run
double FastTimerService::queryModuleTimeByLabel(edm::StreamID sid, const std::string & label) const {
  auto const & keyval = m_stream[sid].modules.find(label);
  if (keyval != m_stream[sid].modules.end()) {
    return keyval->second.time_active;
  } else {
    // module not found
    edm::LogError("FastTimerService") << "FastTimerService::queryModuleTimeByLabel: unexpected module " << label;
    return 0.;
  }
}

// query the time spent in a module (available after the module has run
double FastTimerService::queryModuleTimeByType(edm::StreamID sid, const std::string & type) const {
  auto const & keyval = m_stream[sid].moduletypes.find(type);
  if (keyval != m_stream[sid].moduletypes.end()) {
    return keyval->second.time_active;
  } else {
    // module not found
    edm::LogError("FastTimerService") << "FastTimerService::queryModuleTimeByType: unexpected module type " << type;
    return 0.;
  }
}

// query the time spent in a path (available after the path has run)
double FastTimerService::queryPathActiveTime(edm::StreamID sid, const std::string & path) const {
  PathMap<PathInfo>::const_iterator keyval = m_stream[sid].paths.find(path);
  if (keyval != m_stream[sid].paths.end()) {
    return keyval->second.time_active;
  } else {
    edm::LogError("FastTimerService") << "FastTimerService::queryPathActiveTime: unexpected path " << path;
    return 0.;
  }
}

// query the time spent in a path (available after the path has run)
double FastTimerService::queryPathExclusiveTime(edm::StreamID sid, const std::string & path) const {
  PathMap<PathInfo>::const_iterator keyval = m_stream[sid].paths.find(path);
  if (keyval != m_stream[sid].paths.end()) {
    return keyval->second.time_exclusive;
  } else {
    edm::LogError("FastTimerService") << "FastTimerService::queryPathExclusiveTime: unexpected path " << path;
    return 0.;
  }
}

// query the total time spent in a path (available after the path has run)
double FastTimerService::queryPathTotalTime(edm::StreamID sid, const std::string & path) const {
  PathMap<PathInfo>::const_iterator keyval = m_stream[sid].paths.find(path);
  if (keyval != m_stream[sid].paths.end()) {
    return keyval->second.time_total;
  } else {
    edm::LogError("FastTimerService") << "FastTimerService::queryPathTotalTime: unexpected path " << path;
    return 0.;
  }
}

// query the time spent in the current event's source (available during event processing)
double FastTimerService::querySourceTime(edm::StreamID sid) const {
  return m_event[sid].source;
}

// query the time spent in the current event's paths (available during endpaths)
double FastTimerService::queryPathsTime(edm::StreamID sid) const {
  return m_event[sid].all_paths;
}

// query the time spent in the current event's endpaths (available after all endpaths have run)
double FastTimerService::queryEndPathsTime(edm::StreamID sid) const {
  return m_event[sid].all_endpaths;
}

// query the time spent processing the current event (available after the event has been processed)
double FastTimerService::queryEventTime(edm::StreamID sid) const {
  return m_event[sid].event;
}

// describe the module's configuration
void FastTimerService::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<bool>(   "useRealTimeClock",         true);
  desc.addUntracked<bool>(   "enableTimingPaths",        true);
  desc.addUntracked<bool>(   "enableTimingModules",      true);
  desc.addUntracked<bool>(   "enableTimingExclusive",    false);
  desc.addUntracked<bool>(   "enableTimingSummary",      false);
  desc.addUntracked<bool>(   "skipFirstPath",            false),
  desc.addUntracked<bool>(   "enableDQM",                true);
  desc.addUntracked<bool>(   "enableDQMbyPathActive",    false);
  desc.addUntracked<bool>(   "enableDQMbyPathTotal",     true);
  desc.addUntracked<bool>(   "enableDQMbyPathOverhead",  false);
  desc.addUntracked<bool>(   "enableDQMbyPathDetails",   false);
  desc.addUntracked<bool>(   "enableDQMbyPathCounters",  true);
  desc.addUntracked<bool>(   "enableDQMbyPathExclusive", false);
  desc.addUntracked<bool>(   "enableDQMbyModule",        false);
  desc.addUntracked<bool>(   "enableDQMbyModuleType",    false);
  desc.addUntracked<bool>(   "enableDQMSummary",         false);
  desc.addUntracked<bool>(   "enableDQMbyLuminosity",    false);
  desc.addUntracked<bool>(   "enableDQMbyLumiSection",   false);
  desc.addUntracked<bool>(   "enableDQMbyProcesses",     false);
  desc.addUntracked<double>( "dqmTimeRange",             1000. );   // ms
  desc.addUntracked<double>( "dqmTimeResolution",           5. );   // ms
  desc.addUntracked<double>( "dqmPathTimeRange",          100. );   // ms
  desc.addUntracked<double>( "dqmPathTimeResolution",       0.5);   // ms
  desc.addUntracked<double>( "dqmModuleTimeRange",         40. );   // ms
  desc.addUntracked<double>( "dqmModuleTimeResolution",     0.2);   // ms
  desc.addUntracked<double>( "dqmLuminosityRange",      1.e34  );   // cm-2 s-1
  desc.addUntracked<double>( "dqmLuminosityResolution", 1.e31  );   // cm-2 s-1
  desc.addUntracked<uint32_t>( "dqmLumiSectionsRange",   2500  );   // ~ 16 hours
  desc.addUntracked<std::string>(   "dqmPath",           "HLT/TimerService");
  desc.addUntracked<edm::InputTag>( "luminosityProduct", edm::InputTag("hltScalersRawToDigi"));
  desc.addUntracked<std::vector<unsigned int> >("supportedProcesses", { 8, 24, 32 });
  descriptions.add("FastTimerService", desc);
}
