#include "../interface/LinearizedTrackFitter.h"

LinearizedTrackFitter::LinearizedTrackFitter(const std::string & baseDir, const bool inputExtrapolateR,
                                             const bool inputCorrectNonRadialStrips) :
  preEstimatePtDirName_("PreEstimate_Transverse"),
  preEstimateCotThetaDirName_("PreEstimate_Longitudinal_Rz"),
  linearFitLowPtDirName_("Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10"),
  linearFitHighPtDirName_("Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more"),
  linearFitLongitudinalDirName_("Combinations_Longitudinal_Rz_SecondOrder"),
  preEstimatedPt_(0.),
  ptSplitValue_(10.),
  combinationIndex_(0),
  baseDir_(baseDir),
  extrapolateR_(false),
  correctNonRadialStrips_(false),
  rotationFactor_(0.)
{
  // Read the constants for the full tracker.
  // It needs to read all 9 or 14 regions based on few base names and generate the names for the removed layers.

  // Fill all pre-estimates
  fillMatrices(baseDir_+"/"+preEstimatePtDirName_, "matrixVD_0_pre_chargeOverPt.txt", &chargeOverPtEstimator_);
  // R and z are assumed to have the same number of layers. If not the estimator needs to be modified.
  fillMatrices(baseDir_+"/"+preEstimateCotThetaDirName_, "matrixVD_0_pre_cotTheta.txt", &cotThetaEstimator_);
  fillMatrices(baseDir_+"/"+preEstimateCotThetaDirName_+"_tgTheta", "matrixVD_0_pre_tgTheta.txt", &tgThetaEstimator_);

  // Fill all PCA coefficients for parameters and chi2 estimates
  // std::cout << baseDir_+"/"+linearFitLowPtDirName_ << std::endl;
  // throw;
  fillMatrices(baseDir_+"/"+linearFitLowPtDirName_, "matrixVD_0.txt", &linearFitLowPt_);
  fillMatrices(baseDir_+"/"+linearFitHighPtDirName_, "matrixVD_0.txt", &linearFitHighPt_);
  fillMatrices(baseDir_+"/"+linearFitLongitudinalDirName_, "matrixVD_0.txt", &linearFitLongitudinal_);
}


/// This is used by the full simulation
double LinearizedTrackFitter::fit(const std::vector<double> & vars, const int bits)
{
  std::vector<int> layers;
  // std::cout << "bits = " << bits << std::endl;
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
/*
  if (bits > 0) {
    std::vector<double> cleanedVars;
    for (size_t i = 0; i < vars.size() / 3; ++i) {
      if (int(i) != (bits - 1)) {
        cleanedVars.push_back(vars.at(i * 3));
        cleanedVars.push_back(vars.at(i * 3 + 1));
        cleanedVars.push_back(vars.at(i * 3 + 2));
      }
    }
    // return fit(cleanedVars, layers, 0., 0.);
    return fit(cleanedVars, layers);
  }
*/
  // return fit(vars, layers, 0., 0.);
  return fit(vars, layers);
}


void LinearizedTrackFitter::initialize(const std::vector<double> & vars, const std::vector<int> & layers)
{
  if (vars.size() < 15) {
    std::cout << "Error: number of input variables is less than 15. Please provide 5 or 6 sets of (phi, R, z) ordered from the innermost to the outermost layer." << std::endl;
    std::cout << "Number of input variables = " << vars.size() << std::endl;
    throw;
  }
  if (layers.size()*3 != vars.size()) {
    std::cout << "Error: inconsistent number of layers and number of variables. They should be in a ratio of 1/3." << std::endl;
    std::cout << "Number of layers = " << layers.size() << std::endl;
    for (auto l : layers) std::cout << l << ", ";
    std::cout << std::endl;
    std::cout << "Number of variables = " << vars.size() << std::endl;
  }
  varsNum_ = vars.size()/3;
  varsR_.clear();
  varsR_.reserve(varsNum_);
  correctedVarsPhi_ = Matrix<long double, Dynamic, 1>(varsNum_);
  correctedVarsZ_ = Matrix<long double, Dynamic, 1>(varsNum_);

  computeRotationFactor(vars);

  // Special care is needed for the +/-pi boundary
  for (unsigned int i=0; i<varsNum_; ++i) {
    correctedVarsPhi_(i) = vars[i*3] - rotationFactor_;
    if (fabs(vars[0]) > M_PI_2) {
      if (vars[0] < 0. && vars[i * 3] > 0.) correctedVarsPhi_(i) -= 2 * M_PI;
      else if (vars[0] > 0. && vars[i * 3] < 0.) correctedVarsPhi_(i) += 2 * M_PI;
    }
  }
  for (unsigned int i=0; i<varsNum_; ++i) { varsR_.push_back(vars[i*3+1]); }
  for (unsigned int i=0; i<varsNum_; ++i) { correctedVarsZ_(i) = vars[i*3+2]; }
  extrapolatedR_ = varsR_;

  uniqueLayers_ = layers;
  std::sort(uniqueLayers_.begin(), uniqueLayers_.end());
  uniqueLayers_.erase(std::unique(uniqueLayers_.begin(), uniqueLayers_.end()), uniqueLayers_.end());
  combinationIndex_ = combinationIndex(uniqueLayers_, varsR_);
}


double LinearizedTrackFitter::fit(const std::vector<double> & vars, const std::vector<int> & layers)
{
  initialize(vars, layers);

  auto iterPt = chargeOverPtEstimator_.find(combinationIndex_);
  if (iterPt == chargeOverPtEstimator_.end()) {
    return -1.;
  }

  if (!readMean(baseDir_+"/"+preEstimateCotThetaDirName_, "MeanRadius_", combinationIndex_, meanRadius_)) {
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

  return fit(chargeOverTwoRho, cotTheta, 0.);
}


double LinearizedTrackFitter::fit(const double & chargeOverTwoRho, const double & cotTheta, const double & tgTheta)
{
  for (unsigned int i=0; i<varsNum_; ++i) {
    double DeltaR = varsR_[i] - meanRadius_[combinationIndex_][i];
    double RCube = std::pow(varsR_[i], 3);
    // Note: the extrapolatedR = R and is only extrapolated for 2S modules in the disks if requested
    double DeltaExtrapolatedR = extrapolatedR_[i] - meanRadius_[combinationIndex_][i];
    double extrapolatedRCube = std::pow(extrapolatedR_[i], 3);
    correctedVarsPhi_[i] += chargeOverTwoRho * DeltaExtrapolatedR + extrapolatedRCube * std::pow(chargeOverTwoRho, 3) / 6.;
    // We use the regular R for the R-z plane. We could recompute constants with the extrapolated R (likely negligible difference).
    correctedVarsZ_[i] -= (DeltaR + 1/6.*RCube*(chargeOverTwoRho*chargeOverTwoRho))*cotTheta;
  }

  // Evaluate the chi2/ndof
  MatrixReader * linearFitLongitudinal = &(linearFitLongitudinal_.find(combinationIndex_)->second);
  double normChi2 = linearFitLongitudinal->normChi2(correctedVarsZ_)*linearFitLongitudinal->nDof();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  normChi2 += linearFitTransverse->normChi2(correctedVarsPhi_)*linearFitTransverse->nDof();
  normChi2 /= linearFitLongitudinal->nDof()+linearFitTransverse->nDof();

  // Estimate the track parameters
  estimatedPars_.clear();
  estimatedPars_ = linearFitTransverse->trackParameters(correctedVarsPhi_);
  // Parameter 1 must be phi0 for the rotation. Keep phi0 in the [-pi, pi] range.
  if (estimatedPars_.size() > 1) {
    estimatedPars_.at(1) += rotationFactor_;
    if (estimatedPars_.at(1) > M_PI) estimatedPars_.at(1) -= 2*M_PI;
    else if (estimatedPars_.at(1) < -M_PI) estimatedPars_.at(1) += 2*M_PI;
  }
  auto tempPars = linearFitLongitudinal->trackParameters(correctedVarsZ_);
  estimatedPars_.insert(estimatedPars_.end(), tempPars.begin(), tempPars.end());

  return normChi2;
}


std::vector<double> LinearizedTrackFitter::principalComponents()
{
  principalComponents_.clear();
  MatrixReader * linearFitTransverse = nullptr;
  if (preEstimatedPt_ < ptSplitValue_) linearFitTransverse = &(linearFitLowPt_.find(combinationIndex_)->second);
  else linearFitTransverse = &(linearFitHighPt_.find(combinationIndex_)->second);
  principalComponents_ = linearFitTransverse->principalComponents(correctedVarsPhi_);
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.principalComponents(correctedVarsZ_);
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
  auto tempPrincipalFromZ = linearFitLongitudinal_.find(combinationIndex_)->second.normalizedPrincipalComponents(correctedVarsZ_);
  normalizedPrincipalComponents_.insert(normalizedPrincipalComponents_.end(), tempPrincipalFromZ.begin(), tempPrincipalFromZ.end());
  return normalizedPrincipalComponents_;
}
