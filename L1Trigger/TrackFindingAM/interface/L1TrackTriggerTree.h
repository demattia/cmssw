#ifndef L1TRACKTRIGGERTREE_H
#define L1TRACKTRIGGERTREE_H

#include <vector>
#include <iostream>
#include "TString.h"
#include "TChain.h"
#include <cmath>

class L1TrackTriggerTree
{
public:
  L1TrackTriggerTree(const TString & fileName);

  // Load the values for the event
  void getEntry(const int evt)
  {
    L1TT->GetEntry(evt);
  }

  void printInfo();

  TChain *L1TT;
  int n_entries;
  int m_stub;
  // std::vector<float>  *m_stub_pt;
  // std::vector<float>  *m_stub_ptGEN;
  std::vector<float>  *m_stub_pxGEN;
  std::vector<float>  *m_stub_pyGEN;
  std::vector<float>  *m_stub_etaGEN;
  std::vector<float>  *m_stub_X0;
  std::vector<float>  *m_stub_Y0;
  std::vector<float>  *m_stub_Z0;
  std::vector<float>  *m_stub_PHI0;
  std::vector<int>    *m_stub_layer;
  std::vector<int>    *m_stub_module;
  std::vector<int>    *m_stub_ladder;
  // std::vector<int>    *m_stub_seg;
  // std::vector<int>    *m_stub_modid;
  std::vector<float>    *m_stub_strip;
  std::vector<float>  *m_stub_x;
  std::vector<float>  *m_stub_y;
  std::vector<float>  *m_stub_z;
  // std::vector<int>    *m_stub_clust1;
  // std::vector<int>    *m_stub_clust2;
  // std::vector<int>    *m_stub_cw1;
  // std::vector<int>    *m_stub_cw2;
  std::vector<float>  *m_stub_deltas;
  // std::vector<float>  *m_stub_cor;
  // std::vector<int>    *m_stub_tp;
  std::vector<int>    *m_stub_pdg;
  std::vector<int>    *m_stub_tp;
  // std::vector<int>    *m_stub_pid;

  std::vector<float>  *m_stub_PHI0Extrapolated;
  std::vector<float>  *m_stub_d0GEN;
  std::vector<float>  *m_stub_z0GENExtrapolated;
};

#endif // L1TRACKTRIGGERTREE_H
