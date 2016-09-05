#ifndef GETTRACKPARAMETERS_H
#define GETTRACKPARAMETERS_H

#include <memory>
#include <math.h>
#include <vector>
#include "L1Trigger/TrackFindingAM/interface/L1TrackTriggerTree.h"
#include <random>

// Abstract base class
class GetTreeTrackParameter
{
public:
  virtual ~GetTreeTrackParameter() {}
  virtual double at(const int k) = 0;
private:
};


// Phi parameter of the generated track associated to stub k
class GetParPhiPrompt : public GetTreeTrackParameter
{
public:
  GetParPhiPrompt(std::shared_ptr<L1TrackTriggerTree> tree) : par_phi(tree->m_stub_PHI0) {}
  virtual ~GetParPhiPrompt() {}
  virtual double at(const int k) {return par_phi->at(k);}
private:
  std::vector<float> * par_phi;
};


class GetParOneOverPt : public GetTreeTrackParameter
{
public:
  GetParOneOverPt(std::shared_ptr<L1TrackTriggerTree> tree) : par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN) {}
  virtual ~GetParOneOverPt() {}
  virtual double at(const int k) {
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    return pt > 0 ? 1./pt : 0.;
  }
private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
};


class GetParChargeOverPt : public GetTreeTrackParameter
{
public:
  GetParChargeOverPt(std::shared_ptr<L1TrackTriggerTree> tree) : par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {}
  virtual ~GetParChargeOverPt() {}
  virtual double at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    return (pt > 0 ? charge/pt : 0.);
  }
private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
};


class GetParChargeOverPtEnergyLossCorrected : public GetTreeTrackParameter
{
 public:
  GetParChargeOverPtEnergyLossCorrected(std::shared_ptr<L1TrackTriggerTree> tree) : par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg) {}
  virtual ~GetParChargeOverPtEnergyLossCorrected() {}
  virtual double at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    // Shifting by 4 MeV to compensate for average energy loss
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2)) - 0.004;
    return (pt > 0 ? charge/pt : 0.);
  }
 private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
};


class GetParCharge : public GetTreeTrackParameter
{
public:
  GetParCharge(std::shared_ptr<L1TrackTriggerTree> tree) : par_pdg(tree->m_stub_pdg) {}
  virtual ~GetParCharge() {}
  virtual double at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    return ((par_pdg->at(k) > 0) ? -1 : 1);
  }
private:
  std::vector<int> * par_pdg;
};


// cotTheta parameter of the generated track associated to stub k
class GetParCotTheta : public GetTreeTrackParameter
{
public:
  GetParCotTheta(std::shared_ptr<L1TrackTriggerTree> tree) : par_eta(tree->m_stub_etaGEN) {}
  virtual ~GetParCotTheta() {}
  virtual double at(const int k) {return 1./tan(2*atan(exp(-par_eta->at(k))));}
private:
  std::vector<float> * par_eta;
};


// z0 parameter of the generated track associated to stub k
class GetParZ0 : public GetTreeTrackParameter
{
public:
  GetParZ0(std::shared_ptr<L1TrackTriggerTree> tree) : par_z0(tree->m_stub_z0GENExtrapolated) {}
  virtual ~GetParZ0() {}
  virtual double at(const int k) {return par_z0->at(k);}
private:
  std::vector<float> * par_z0;
};


// z0 parameter of the generated track associated to stub k
class GetParZ0Prompt : public GetTreeTrackParameter
{
public:
  GetParZ0Prompt(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_z0(tree->m_stub_Z0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0), par_eta(tree->m_stub_etaGEN)  {}
  virtual ~GetParZ0Prompt() {}
  virtual double at(const int k) {
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    double r = pt / (0.003 * 3.8114); // In centimeters (0.3 for meters)
    double xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    double yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    // The impact parameter is the distance between the trajectory (simplified as a circle in the transverse plane)
    // and the origin (which is the reference point in this case). It will need to be adapted for a beamspot.
    double d = r - std::sqrt(xc*xc + yc*yc);
    double cotTheta = 1./tan(2*atan(exp(-par_eta->at(k))));
    return par_z0->at(k) - (std::sqrt(par_x0->at(k)*par_x0->at(k) + par_y0->at(k)*par_y0->at(k)) - d)*cotTheta;
  }
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_z0;
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_phi;
  std::vector<float> * par_eta;
};


// d0 parameter of the generated track associated to stub k
class GetParD0 : public GetTreeTrackParameter
{
public:
  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0), par_d0(tree->m_stub_d0GEN) {}
  virtual ~GetParD0() {}
  virtual double at(const int k) {
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
    double r = pt / (0.003 * 3.8114); // In centimeters (0.3 for meters)
    double xc = -charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    double yc = charge*r * cos(par_phi->at(k)) + par_y0->at(k);
    // The impact parameter is the distance between the trajectory (simplified as a circle in the transverse plane)
    // and the origin (which is the reference point in this case). It will need to be adapted for a beamspot.
    // The charge is needed so that it reflects consistently with the variables.
    double d0Sign = charge*(r - std::sqrt(xc*xc + yc*yc));
    return (d0Sign < 0 ? -par_d0->at(k) : par_d0->at(k));
  }
private:
  std::vector<float> * par_x0;
  std::vector<float> * par_y0;
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_phi;
  std::vector<float> * par_d0;
};


// phi0 parameter of the generated track associated to stub k at the point of closest approach
class GetParPhi : public GetTreeTrackParameter
{
 public:
  GetParPhi(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_phi0(tree->m_stub_PHI0Extrapolated) {}
  virtual ~GetParPhi() {}
  virtual double at(const int k) {
//    std::cout << "extrapolated phi0 = " << par_phi0->at(k) << std::endl;
    return par_phi0->at(k);
  }
 private:
  std::vector<float> * par_phi0;
};


// Endcap parameters
class GetParZ0TgTheta : public GetTreeTrackParameter
{
 public:
  GetParZ0TgTheta(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_z0(tree->m_stub_Z0), par_eta(tree->m_stub_etaGEN) {}
  virtual ~GetParZ0TgTheta() {}
  virtual double at(const int k) {
    return par_z0->at(k)*tan(2*atan(exp(-par_eta->at(k))));
  }
 private:
  std::vector<float> * par_z0;
  std::vector<float> * par_eta;
};


class GetParTgTheta : public GetTreeTrackParameter
{
 public:
  GetParTgTheta(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_eta(tree->m_stub_etaGEN) {}
  virtual ~GetParTgTheta() {}
  virtual double at(const int k) {
    return tan(2*atan(exp(-par_eta->at(k))));
  }
 private:
  std::vector<float> * par_eta;
};


class GetParChargeOverPz : public GetTreeTrackParameter
{
 public:
  GetParChargeOverPz(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg), par_eta(tree->m_stub_etaGEN) {}
  virtual ~GetParChargeOverPz() {}
  virtual double at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pz = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2))/tan(2*atan(exp(-par_eta->at(k))));
    return (pz > 0 ? charge/pz : 0.);
  }
 private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_eta;
};


class GetParPhi0PlusChargeZ0Over2RhoZ : public GetTreeTrackParameter
{
 public:
  GetParPhi0PlusChargeZ0Over2RhoZ(std::shared_ptr<L1TrackTriggerTree> tree) :
      par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN), par_pdg(tree->m_stub_pdg),
      par_z0(tree->m_stub_Z0), par_eta(tree->m_stub_etaGEN), par_phi0(tree->m_stub_PHI0Extrapolated) {}
  virtual ~GetParPhi0PlusChargeZ0Over2RhoZ() {}
  virtual double at(const int k) {
    // For muons, electrons and taus the charge is the opposite of the sign of the pdgId
    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
    double pz = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2))/tan(2*atan(exp(-par_eta->at(k))));
    return (pz > 0 ? par_phi0->at(k) + par_z0->at(k)*3.8114*0.003/2.*charge/pz : 0.);
  }
 private:
  std::vector<float> * par_px;
  std::vector<float> * par_py;
  std::vector<int> * par_pdg;
  std::vector<float> * par_z0;
  std::vector<float> * par_eta;
  std::vector<float> * par_phi0;
};

#endif // GETTRACKPARAMETERS_H
