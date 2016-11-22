#include "L1Trigger/TrackFindingAM/interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                                             const int inputExtrapolatedRPrecision,
                                             const bool inputCorrectNonRadialStrips, const int regionsNumber,
                                             const bool doCutOnPrincipals,
                                             const std::string & preEstimatePtDirName,
                                             const std::string & preEstimateCotThetaDirName,
                                             const std::string & linearFitLowPtDirName,
                                             const std::string & linearFitHighPtDirName,
                                             const std::string & linearFitLongitudinalDirName,
                                             const bool alignPrincipals) :
    LinearizedTrackFitterBase(baseDir, inputExtrapolateR, inputExtrapolatedRPrecision,
                              inputCorrectNonRadialStrips, regionsNumber, doCutOnPrincipals,
                              preEstimatePtDirName, preEstimateCotThetaDirName,
                              linearFitLowPtDirName, linearFitHighPtDirName, linearFitLongitudinalDirName,
                              alignPrincipals),
    preEstimatedPt_(0.),
    cuts6Low_{-1., 45., 35., 25., -1., -1., -1, 3., 3., 3., -1., -1.},
    cuts6High_{-1., 55., 45., 40., -1., -1., -1, 3., 3., 3., -1., -1.}
{
  // Fill all pre-estimates
  fillMatrices(preEstimatePtDirName_, "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_);
  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_);
  fillMatrices(preEstimateCotThetaDirName_, "matrixVD_0_pre_tgTheta.txt", &tgThetaEstimator_);

  // Fill all PCA coefficients for parameters and chi2 estimates
  fillMatrices(linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_);
  fillMatrices(linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_);
  fillMatrices(linearFitLongitudinalDirName_, "matrixVD_0.txt", &linearFitLongitudinal_);
}


/// This is used by the full simulation
double LinearizedTrackFitter::fit(const std::vector<double> & vars, const int bits)
{
  std::vector<int> layers;
  if (bits == 0) layers = {5, 6, 7, 8, 9, 10};
  else if (bits == 1) layers = {6, 7, 8, 9, 10};
  else if (bits == 2) layers = {5, 7, 8, 9, 10};
  else if (bits == 3) layers = {5, 6, 8, 9, 10};
  else if (bits == 4) layers = {5, 6, 7, 9, 10};
  else if (bits == 5) layers = {5, 6, 7, 8, 10};
  else if (bits == 6) layers = {5, 6, 7, 8, 9};
  else {
    std::cout << "Error: unknown bits = " << bits << std::endl;
    throw;
  }

  // Clean the variables removing the 0 values corresponding to the missing layer
  if (bits > 0) {
    std::vector<double> cleanedVars;
    for (size_t i = 0; i < vars.size() / 3; ++i) {
      if (int(i) != (bits - 1)) {
        cleanedVars.push_back(vars.at(i * 3));
        cleanedVars.push_back(vars.at(i * 3 + 1));
        cleanedVars.push_back(vars.at(i * 3 + 2));
      }
    }
    return fit(cleanedVars, layers);
  }

  return fit(vars, layers);
}


bool LinearizedTrackFitter::initialize(const std::vector<double> & vars, const std::vector<int> & layers)
{
  bool validInput = consistencyCheck(vars, layers);
  if (!validInput) return false;
  varsNum_ = vars.size()/3;
  varsR_.clear();
  varsR_.reserve(varsNum_);
  correctedVarsPhi_ = Matrix<long double, Dynamic, 1>(varsNum_);
  correctedVarsZ_ = Matrix<long double, Dynamic, 1>(varsNum_);

  // Compute the phi alignment value
  rotationFactor_ = computeRotationFactor(vars);

  // If there are stubs in the negative disks flip the z coordinates (use z symmetry)
  // negativeZ_ = false;
  // for (auto l : layers) {
  //   if (l > 15) {
  //     negativeZ_ = true;
  //   }
  // }

  // Invert the z coordinates if all the z are negative. If any of the z is positive do nothing.
  negativeZ_ = true;
  for (unsigned int i=0; i<varsNum_; ++i) {
    if (vars[i*3+2] > 0.) negativeZ_ = false;
  }

  // Rotate phi accounting for the discontinuity at +/-pi
  for (unsigned int i=0; i<varsNum_; ++i) {
    if (rotationFactor_ < -1.2 && vars[i*3] > 0.) correctedVarsPhi_(i) = vars[i*3]-2*M_PI - rotationFactor_;
    else if (rotationFactor_ > 1.2 && vars[i*3] < 0.) correctedVarsPhi_(i) = vars[i*3]+2*M_PI - rotationFactor_;
    else correctedVarsPhi_(i) = vars[i*3] - rotationFactor_;
  }
  for (unsigned int i=0; i<varsNum_; ++i) { varsR_.push_back(vars[i*3+1]); }
  for (unsigned int i=0; i<varsNum_; ++i) { correctedVarsZ_(i) = (negativeZ_ ? -vars[i*3+2] : vars[i*3+2]); }
  extrapolatedR_ = varsR_;

  uniqueLayers_ = layers;
  // Transform all negative disks into the corresponding positive disks.
  if (negativeZ_) {
    std::transform(uniqueLayers_.begin(), uniqueLayers_.end(), uniqueLayers_.begin(), [](int l) -> int { return (l > 15 ? l - 5 : l); });
    // std::cout << "inverted z layers = ";
    // for (auto l : uniqueLayers_) std::cout << l << " ";
    // std::cout << std::endl;
  }
  std::sort(uniqueLayers_.begin(), uniqueLayers_.end());
  uniqueLayers_.erase(std::unique(uniqueLayers_.begin(), uniqueLayers_.end()), uniqueLayers_.end());
  combinationIndex_ = combinationIndex(uniqueLayers_, varsR_, regionsNumber_);
  return true;
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars, const std::vector<int> & layers)
{
  // std::cout << "layers: " << std::endl;
  // for (auto l : layers) std::cout << l << ", " << std::endl;
  // std::cout << std::endl;
  bool validInput = initialize(vars, layers);
  if (!validInput) {
    std::cout << "invalid input" << std::endl;
    // throw;
    return -1.;
  }

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    // std::cout << "missing combination" << std::endl;
    return -1.;
  }

  if (!readMean(preEstimateCotThetaDirName_, "MeanRadius_", combinationIndex_, meanRadius_)) {
    std::cout << "Error: mean radii not found for combination = " << combinationIndex_ << std::endl;
    throw;
  }

  // Correct the input variables and split them between phi and z vectors
  EstimatorSimple & chargeOverPtEstimator = iterPt->second;
  double preEstimatedChargeOverPt = chargeOverPtEstimator.estimate(correctedVarsPhi_);
  auto iterCotTheta = cotThetaEstimator_.find(combinationIndex_);
  EstimatorSimple & cotThetaEstimator = iterCotTheta->second;

  preEstimatedPt_ = 1./fabs(preEstimatedChargeOverPt);
  // Retake it here because we need it with the charge
  double chargeOverTwoRho = (3.8114*0.003)*preEstimatedChargeOverPt/2.;
  double cotTheta = cotThetaEstimator.estimate(varsR_, correctedVarsZ_);

  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  double tgTheta = 0.;
  if (extrapolateR_) {
    auto iterTgTheta = tgThetaEstimator_.find(combinationIndex_);
    EstimatorSimple &tgThetaEstimator = iterTgTheta->second;
    tgTheta = tgThetaEstimator.estimate(varsR_, correctedVarsZ_);
  }

  return fit(chargeOverTwoRho, cotTheta, tgTheta);
}


double LinearizedTrackFitter::fit(const double & chargeOverTwoRho, const double & cotTheta, const double & tgTheta)
{
  // Extrapolate R if required
  // Warning: do not put in the following loop or the correctedVarsZ will be modified for all the elements after the
  // first one and the results will be incorrect.
  if (extrapolateR_) {
    firstOrderTerm_.clear();
    secondOrderTerm1_.clear();
    secondOrderTerm2_.clear();
    secondOrderTerm3_.clear();

    // Force first order R extrapolation
    // extrapolatedRPrecision_ = 0;

    for (unsigned int i=0; i<varsNum_; ++i) {
      if (extrapolatedRPrecision_ == 0) {
        extrapolatedR_[i] = extrapolateRFirstOrder(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                   chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                   firstOrderTerm_);
      }
      else if (extrapolatedRPrecision_ == 1) {
        extrapolatedR_[i] = extrapolateRSecondOrderFirstTermOnly(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                                 chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                                 firstOrderTerm_, secondOrderTerm1_);
      }
      else if (extrapolatedRPrecision_ == 2) {
        extrapolatedR_[i] = extrapolateRSecondOrderFirstTwoTermsOnly(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                                     chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                                     firstOrderTerm_, secondOrderTerm1_, secondOrderTerm2_);
      }
      else if (extrapolatedRPrecision_ == 3) {
        extrapolatedR_[i] = extrapolateRSecondOrder(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                                    chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_,
                                                    firstOrderTerm_, secondOrderTerm1_, secondOrderTerm2_,
                                                    secondOrderTerm3_);
      }
      else if (extrapolatedRPrecision_ == 4) {
        extrapolatedR_[i] = extrapolateRExact(varsR_[i], correctedVarsZ_[i], uniqueLayers_[i], tgTheta,
                                              chargeOverTwoRho, uniqueLayers_, varsR_, correctedVarsZ_);
      }
      if (correctNonRadialStrips_) {
        // We correct for the rotation factor to allow for the lookup table to work properly.
        correctedVarsPhi_[i] = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(correctedVarsPhi_[i]+rotationFactor_, 0.009, extrapolatedR_[i],
                                                                                                varsR_[i], correctedVarsZ_[i], uniqueLayers_[i]) - rotationFactor_;
      }
    }
  }

  for (unsigned int i=0; i<varsNum_; ++i) {
    double DeltaR = varsR_[i] - meanRadius_[combinationIndex_][i];
    double RCube = std::pow(varsR_[i], 3);
    // Note: the extrapolatedR = R and is only extrapolated for 2S modules in the disks if requested
    double DeltaExtrapolatedR = extrapolatedR_[i] - meanRadius_[combinationIndex_][i];
    double extrapolatedRCube = std::pow(extrapolatedR_[i], 3);
    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaExtrapolatedR + extrapolatedRCube * std::pow(chargeOverTwoRho, 3) / 6.;
//    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaExtrapolatedR;
//    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaR + RCube * std::pow(chargeOverTwoRho, 3) / 6.;
//    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaR;
    // We use the regular R for the R-z plane. We could recompute constants with the extrapolated R (likely negligible difference).
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(chargeOverTwoRho*chargeOverTwoRho))*cotTheta;
//    correctedVarsZ_[i] -= DeltaR * cotTheta;
    // correctedVarsZ_[i] -= DeltaExtrapolatedR*cotTheta;
  }

  // Evaluate the chi2/ndof
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  ndofLongitudinal_ = linearFitLongitudinal->nDof();
  chi2Longitudinal_ = linearFitLongitudinal->normChi2(correctedVarsZ_)*ndofLongitudinal_;
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  ndofTransverse_ = linearFitTransverse->nDof();
  chi2Transverse_ = linearFitTransverse->normChi2(correctedVarsPhi_)*ndofTransverse_;

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhi_);
  // Parameter 1 must be phi0 for the rotation.
  if (estimatedPars_.size() > 1) estimatedPars_.at(1) += rotationFactor_;
  auto tempPars = linearFitLongitudinal->trackParameters(correctedVarsZ_);
  // Change the sign of cot(theta) and z0 if we originally changed the sign of the z coordinates
  // Note that this requires for cot(theta) and z0 to be the first two parameters returned by the longitudinal fit.
  if (negativeZ_ && (tempPars.size() > 1)) {
    tempPars.at(0) = -tempPars.at(0);
    tempPars.at(1) = -tempPars.at(1);
  }
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());

  if (doCutOnPrincipals_) {
    if (!principalCuts(normalizedPrincipalComponents(), preEstimatedPt_)) return -1;
  }

  return (chi2Transverse_+chi2Longitudinal_)/(ndofTransverse_+ndofLongitudinal_);
}


std::vector<double> LinearizedTrackFitter::principalComponents()
{
  principalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  if (alignPrincipals_) alignPrincipals(principalComponents_, linearFitTransverse->nDof());
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  auto tempPrincipalFromZ = linearFitLongitudinal->principalComponents(correctedVarsZ_);
  if (alignPrincipals_) alignPrincipals(tempPrincipalFromZ, linearFitLongitudinal->nDof());
  principalComponents_.insert(principalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return principalComponents_;
}


std::vector<double> LinearizedTrackFitter::normalizedPrincipalComponents()
{
  normalizedPrincipalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  normalizedPrincipalComponents_ = linearFitTransverse->normalizedPrincipalComponents(correctedVarsPhi_);
  if (alignPrincipals_) alignPrincipals(normalizedPrincipalComponents_, linearFitTransverse->nDof());
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  auto tempPrincipalFromZ = linearFitLongitudinal->normalizedPrincipalComponents(correctedVarsZ_);
  if (alignPrincipals_) alignPrincipals(tempPrincipalFromZ, linearFitLongitudinal->nDof());
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}


// principal component cut function
bool LinearizedTrackFitter::principalCuts(const std::vector<double> & princes, const float & preEstimatePt)
{
//  std::cout << "preEstimatePt = " << preEstimatePt << std::endl;
//  std::cout << "principals: ";
//  for (auto p : princes) {
//    std::cout << p << ", ";
//  }
//  std::cout << std::endl;
  std::vector<float> * cuts6 = nullptr;
  // Switch cuts based on preEstimatePt
  if (preEstimatePt < 10.) cuts6 = &cuts6Low_;
  else cuts6 = &cuts6High_;
  if (princes.size() == 12) {
    for (unsigned i = 0; i < princes.size(); ++i) {
      if (cuts6->at(i) == -1) continue;
      if (fabs(princes[i]) > cuts6->at(i)) {
//        std::cout << "Fails cut on component " << i << ":" << std::endl;
//        std::cout << "fabs(princes[i]) > cuts6->at(i) = " << fabs(princes[i]) << " > " << cuts6->at(i) << std::endl;
        return false;
      }
    }
  }
  else {
    std::cout << "Error: not the right format!!!!" << std::endl;
    throw;
  }
  return true;
}


//bool LinearizedTrackFitter::principalCuts(const std::vector<double> & princes)
//{
//  bool pass=true;
//  const float Cuts6[12]={58.,33.,22.,12.,-1.,-1.,9.,3.,3.,3.,-1.,-1.};
//  const float Cuts5[10]={33.,22.,12.,-1.,-1.,3.,3.,3.,-1.,-1.};
//  if (princes.size()==12) {
//    for (unsigned i = 0; i < princes.size(); ++i) {
//      if (Cuts6[i] == -1) continue;
//      if (fabs(princes[i]) > Cuts6[i]) {
//        pass = false;
//        break;
//      }
//    }
//  }
//  else {
//    for (unsigned i = 0; i < princes.size(); ++i) {
//      if (Cuts5[i] == -1) continue;
//      if (fabs(princes[i]) > Cuts5[i]) {
//        pass = false;
//        break;
//      }
//    }
//  }
//  return pass;
//}
