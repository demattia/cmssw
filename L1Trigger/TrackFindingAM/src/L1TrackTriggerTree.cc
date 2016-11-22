#include "L1Trigger/TrackFindingAM/interface/L1TrackTriggerTree.h"

L1TrackTriggerTree::L1TrackTriggerTree(const TString & fileName)
{
  L1TT = new TChain("TkStubs");
  L1TT->Add(fileName);

  // L1TT->SetBranchAddress("L1TkSTUB_clust1",    &m_stub_clust1);
  // L1TT->SetBranchAddress("L1TkSTUB_clust2",    &m_stub_clust2);
  // L1TT->SetBranchAddress("L1TkSTUB_cor",       &m_stub_cor);
  L1TT->SetBranchAddress("L1TkSTUB_PHI0",      &m_stub_PHI0);
  // L1TT->SetBranchAddress("L1TkSTUB_tp",        &m_stub_tp);
  L1TT->SetBranchAddress("L1TkSTUB_pdgID",     &m_stub_pdg);
  L1TT->SetBranchAddress("L1TkSTUB_tp",        &m_stub_tp);
  // L1TT->SetBranchAddress("L1TkSTUB_process",   &m_stub_pid);
  L1TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  // L1TT->SetBranchAddress("L1TkSTUB_pt",        &m_stub_pt);
  // L1TT->SetBranchAddress("L1TkSTUB_ptGEN",     &m_stub_ptGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pxGEN",     &m_stub_pxGEN);
  L1TT->SetBranchAddress("L1TkSTUB_pyGEN",     &m_stub_pyGEN);
  L1TT->SetBranchAddress("L1TkSTUB_etaGEN",    &m_stub_etaGEN);
  L1TT->SetBranchAddress("L1TkSTUB_layer",     &m_stub_layer);
  L1TT->SetBranchAddress("L1TkSTUB_module",    &m_stub_module);
  L1TT->SetBranchAddress("L1TkSTUB_ladder",    &m_stub_ladder);
  // L1TT->SetBranchAddress("L1TkSTUB_seg",       &m_stub_seg);
  // L1TT->SetBranchAddress("L1TkSTUB_modid",     &m_stub_modid);
  L1TT->SetBranchAddress("L1TkSTUB_strip",     &m_stub_strip);
  L1TT->SetBranchAddress("L1TkSTUB_x",         &m_stub_x);
  L1TT->SetBranchAddress("L1TkSTUB_y",         &m_stub_y);
  L1TT->SetBranchAddress("L1TkSTUB_z",         &m_stub_z);
  L1TT->SetBranchAddress("L1TkSTUB_deltas",    &m_stub_deltas);
  L1TT->SetBranchAddress("L1TkSTUB_X0",        &m_stub_X0);
  L1TT->SetBranchAddress("L1TkSTUB_Y0",        &m_stub_Y0);
  L1TT->SetBranchAddress("L1TkSTUB_Z0",        &m_stub_Z0);

  L1TT->SetBranchAddress("L1TkSTUB_PHI0Extrapolated",  &m_stub_PHI0Extrapolated);
  L1TT->SetBranchAddress("L1TkSTUB_d0GEN",             &m_stub_d0GEN);
  L1TT->SetBranchAddress("L1TkSTUB_z0GENExtrapolated", &m_stub_z0GENExtrapolated);

  n_entries = L1TT->GetEntries();
}


void L1TrackTriggerTree::printInfo()
{
  std::cout << "where " << m_stub << " stub(s) where produced" << std::endl;
}
