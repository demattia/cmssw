#ifndef _PATTERNFINDER_H_
#define _PATTERNFINDER_H_

#include <iostream>
#include <fstream>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <memory>
#include "SectorTree.h"
#include "LocalToGlobalConverter.h"
#include "LinearizedTrackFitter.h"

#ifdef IPNL_USE_CUDA
#include "gpu.h"
#endif

using namespace std;

/**
   \brief Used to find active patterns in events. Stubs of the event are injected in a virtual detector, then we count the number of active patterns in the bank (the patterns have previously been linked to the virtual detector's supertrips).
**/
class PatternFinder{
 private:
  int active_threshold;
  int max_nb_missing_hit;
  bool useMissingHits;
  unsigned int max_road_number;
  SectorTree* sectors;
  string eventsFilename;
  string outputFileName;
  Detector tracker;
  LocalToGlobalConverter* converter;

#ifdef IPNL_USE_CUDA
  deviceDetector* d_detector;  
  patternBank* d_p_bank;
  deviceParameters* d_parameters;
  int nb_blocks;
  int nb_threads;
#endif

 public:
 /**
     \brief Constructor
     \param at The minimum number of hit super strip to activate a pattern
     \param st The SectorTree containing the sectors with their associated patterns
     \param f The name of the file to analyse
     \param of The name of the output file
  **/
  PatternFinder(int at, SectorTree* st, string f, string of);

  /**
     \brief Destructor
  **/
  ~PatternFinder();

#ifdef IPNL_USE_CUDA
 /**
     \brief Constructor
     \param at The minimum number of hit super strip to activate a pattern
     \param st The SectorTree containing the sectors with their associated patterns
     \param f The name of the file to analyse
     \param of The name of the output file
     \param p The device pattern bank
     \param d The device detector
     \param d_p Structure containing device addresses where parameters are stored
  **/
  PatternFinder(int at, SectorTree* st, string f, string of, patternBank* p, deviceDetector* d, deviceParameters* dp);


  /**
     \brief Get active patterns from list of hits (public for CMSSW).
  **/
  int findCuda(int nb,  deviceStubs* d_stubs, cudaStream_t* stream=NULL);

#endif

  /**
     \brief Set the SectorTree (contains sectors with their patterns)
     \param s The SectorTree containing the sectors with their associated patterns
  **/
  void setSectorTree(SectorTree* s);

  /**
     \brief Set the maximum output road number
     \param m The number of active roads will be limited to the first m patterns (ordered by popularity)
  **/
  void setMaxRoadNumber(unsigned int m);

  /**
     \brief Set the name of the root file containing events
     \param f The name of the file
  **/
  void setEventsFile(string f);
  /**
     \brief Look for active patterns in events
     \param start The search will start from this event number
     \param stop The search will end at this event number
  **/
  void find(int start, int& stop);

#ifdef IPNL_USE_CUDA
  /**
     \brief Look for active patterns in events
     \param start The search will start from this event number
     \param stop The search will end at this event number
  **/
  void findCuda(int start, int& stop, deviceStubs* d_stubs);
#endif

  /**
     \brief Get active patterns from list of hits (public for CMSSW).
  **/
  vector<Sector*> find(vector<Hit*> hits);

  /**
     \brief Control the output level
     \param m If set to True, all stub's superstrip values will be displayed during the pattern recognition process.
   **/
  void setVerboseMode(bool m);
  /**
     \brief Use the maximum missing hit threshold instead of the active_threshold
   **/
  void useMissingHitThreshold(int max_nb_missing_hit);
};
#endif
