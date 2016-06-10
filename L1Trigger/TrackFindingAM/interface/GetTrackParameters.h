#ifndef GETTRACKPARAMETERS_H
#define GETTRACKPARAMETERS_H

#include <memory>
#include <math.h>
#include <vector>
#include "L1TrackTriggerTree.h"
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
  GetParZ0(std::shared_ptr<L1TrackTriggerTree> tree) : par_z0(tree->m_stub_Z0) {}
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


//// d0 parameter of the generated track associated to stub k
//class GetParD0 : public GetTreeTrackParameter
//{
//public:
//  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) : par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0) {}
//  virtual ~GetParD0() {}
//  virtual double at(const int k) {return std::sqrt(std::pow(par_x0->at(k), 2) + std::pow(par_y0->at(k), 2));}
//private:
//  std::vector<double> * par_x0;
//  std::vector<double> * par_y0;
//};


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
//    double xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
//    double yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);

    double xc = -charge*r * sin(par_phi->at(k)) + par_x0->at(k);
    double yc = charge*r * cos(par_phi->at(k)) + par_y0->at(k);

//    double xc = par_x0->at(k) - r*sin(par_phi->at(k));
//    double yc = par_y0->at(k) + charge*r*cos(par_phi->at(k));

    // The impact parameter is the distance between the trajectory (simplified as a circle in the transverse plane)
    // and the origin (which is the reference point in this case). It will need to be adapted for a beamspot.
//    std::cout << "xc = " << xc << ", yc = " << yc << std::endl;
//    std::cout << "r = " << r << ", d0 = " << fabs(r - std::sqrt(xc*xc + yc*yc)) << std::endl;
//    return fabs(r - std::sqrt(xc*xc + yc*yc));
//    return std::sqrt(xc*xc + yc*yc) - r;
    // The charge is needed so that it reflects consistently with the variables.
//    std::cout << "genD0 = " << par_d0->at(k) << std::endl;
//    std::cout << "d0 = " << charge*(std::sqrt(xc*xc + yc*yc) - r) << std::endl;

      double d0Sign = charge*(r - std::sqrt(xc*xc + yc*yc));
//    double d0Sign = charge*(std::sqrt(xc*xc + yc*yc) - r);
//    return d0Sign;
//    double d0Sign = (std::sqrt(xc*xc + yc*yc) - r);

//    std::cout << "par_x0->at("<<k<<") = " << par_x0->at(k) << std::endl;
//    // std::cout << "par_y0->at("<<k<<") = " << par_y0->at(k) << std::endl;
//    // std::cout << "par_d0->at("<<k<<") = " << par_d0->at(k) << std::endl;
//    std::cout << "signed par_d0->at("<<k<<") = " << (d0Sign < 0 ? -par_d0->at(k) : par_d0->at(k)) << std::endl;


//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> dis(-1, 1);
//    return dis(gen);


    return (d0Sign < 0 ? -par_d0->at(k) : par_d0->at(k));


    // return sqrt(par_x0->at(k)*par_x0->at(k) + par_y0->at(k)*par_y0->at(k));

    //return (-par_x0->at(k)*sin(par_phi->at(k))+par_y0->at(k)*cos(par_phi->at(k)));


//    return par_d0->at(k);
    // double Rc = std::sqrt(xc*xc + yc*yc);
    // double xd = xc*(1-r/Rc);
    // double yd = yc*(1-r/Rc);
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


//// d0 parameter of the generated track associated to stub k
//class GetParD0 : public GetTreeTrackParameter
//{
// public:
//  GetParD0(std::shared_ptr<L1TrackTriggerTree> tree) :
//      // par_pdg(tree->m_stub_pdg),
//      par_d0(tree->m_stub_d0GEN) {}
//  virtual ~GetParD0() {}
//  virtual double at(const int k) {
//    // int charge = (par_pdg->at(k) > 0 ? -1 : 1);
//    // return charge*par_d0->at(k);
//    return par_d0->at(k);
//  }
// private:
//  // std::vector<int> * par_pdg;
//  std::vector<double> * par_d0;
//};


////// phi0 parameter of the generated track associated to stub k at the point of closest approach
//class GetParPhi : public GetTreeTrackParameter
//{
//public:
//  GetParPhi(std::shared_ptr<L1TrackTriggerTree> tree) :
//      par_x0(tree->m_stub_X0), par_y0(tree->m_stub_Y0), par_px(tree->m_stub_pxGEN), par_py(tree->m_stub_pyGEN),
//      par_pdg(tree->m_stub_pdg), par_phi(tree->m_stub_PHI0), par_phi0(tree->m_stub_PHI0Extrapolated)
//  {}
//  virtual ~GetParPhi() {}
//  virtual double at(const int k) {
//    int charge = (par_pdg->at(k) > 0 ? -1 : 1);
//    double pt = std::sqrt(std::pow(par_px->at(k), 2) + std::pow(par_py->at(k), 2));
//    double r = pt / (0.003 * 3.8114); // In centimeters (0.3 for meters)
////    double xc = charge*r * sin(par_phi->at(k)) + par_x0->at(k);
////    double yc = -charge*r * cos(par_phi->at(k)) + par_y0->at(k);
////    std::cout << "pt = " << pt << std::endl;
////    std::cout << "vx = " << par_x0->at(k) << ", vy = " << par_y0->at(k) << std::endl;
////    std::cout << "phi0 = " << par_phi->at(k) << std::endl;
//    double xc = -charge*r * sin(par_phi->at(k)) + par_x0->at(k);
//    double yc = charge*r * cos(par_phi->at(k)) + par_y0->at(k);
//
////    double xc = par_x0->at(k) - r*sin(par_phi->at(k));
////    double yc = par_y0->at(k) + charge*r*cos(par_phi->at(k));
//
//
//    // Do not need this because we are computing everything with respect to the origin.
//    double Rc = std::sqrt(xc*xc + yc*yc);
//    double xd = xc*(1-r/Rc);
//    double yd = yc*(1-r/Rc);
//    // double phi_corr = std::atan2(xc-xd, -(yc-yd));
////    if (yc-yd == 0.) return (xc-xd >
//    // To compute the correct angle we need to account for the full sign of the charge.
//    // The angle is defined by the px and py, when using xc and yc to compute it we
//    // need to correct for the difference in sign.
////    double phi_corr = std::atan2(charge*(xc-xd), -charge*(yc-yd));
////    double phi_corr = std::atan2(xc-xd, -(yc-yd));
////    double phi_corr = std::atan2(-charge*(xc-xd), charge*(yc-yd));
////    std::cout << "alpha = " << M_PI_2 - atan2(yc, xc) << std::endl;
////    std::cout << "phi_0' = " << par_phi->at(k) - (M_PI_2 - atan2(yc, xc)) << std::endl;
////    std::cout << "phi_0' true = " << atan2(yc, xc) + M_PI_2 << std::endl;
//    // std::cout << "Extrapolated phi0 = " << par_phi0->at(k) << std::endl;
//    double phi_corr = std::atan2(-charge*(xc-xd), charge*(yc-yd));
//    // std::cout << "phi corr = " << phi_corr << std::endl;
//
////    double phi_corr = std::atan2(-yc, -xc);
//
////    std::cout << "recomputed phi0 = " << atan2(charge*r * cos(par_phi->at(k)), -charge*r * sin(par_phi->at(k)))+M_PI_2 << std::endl;
////    double phi0_recomputed = atan2(charge*r * cos(par_phi->at(k)), charge*r * sin(par_phi->at(k)))+M_PI_2;
////    if (phi0_recomputed < M_PI) phi0_recomputed
////    std::cout << "recomputed phi0 more = " <<  << std::endl;
//
////    double phi_corr = -std::atan(xc/yc);
////    if (phi_corr < 0) phi_corr = phi_corr + M_PI;
////    if (phi_corr > M_PI) phi_corr = phi_corr - M_PI;
////    else if (phi_corr < -M_PI) phi_corr = phi_corr + M_PI;
//    // The minus sign comes from the definition of the angle.
////    std::cout << "phi = " << par_phi->at(k) << std::endl;
////    std::cout << "charge = " << charge << std::endl;
////    std::cout << "x0 = " << par_x0->at(k) << ", y0 = " << par_x0->at(k) << std::endl;
////    std::cout << "xc = " << xc << ", yc = " << yc << std::endl;
//////    std::cout << "xd = " << xd << ", yd = " << yd << std::endl;
////    // std::cout << "phi corr = " << -std::atan2(xc-xd, yc-yd) << std::endl;
////    std::cout << "phi corr = " << phi_corr << std::endl;
//////    return (phi_corr < -M_PI ? phi_corr + M_PI : phi_corr);
////    std::cout << "d0-corrected phi0 = " << phi_corr << std::endl;
//    return phi_corr;
//  }
//private:
//  std::vector<double> * par_x0;
//  std::vector<double> * par_y0;
//  std::vector<double> * par_px;
//  std::vector<double> * par_py;
//  std::vector<int> * par_pdg;
//  std::vector<double> * par_phi;
//  std::vector<double> * par_phi0;
//};


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
