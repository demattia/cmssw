//
// Created by Marco De Mattia on 9/7/15.
//

#ifndef REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H
#define REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H

#include <string>
#include <vector>
#include "L1Trigger/TrackFindingAM/interface/GetVariables.h"
#include "L1Trigger/TrackFindingAM/interface/MatrixReader.h"
#include "L1Trigger/TrackFindingAM/interface/CombinationIndexListBuilder.h"
#include "L1Trigger/TrackFindingAM/interface/BuildTestFunctions.h"

class LinearizedTrackFitterBase
{
 public:
  LinearizedTrackFitterBase(const std::string & baseDir, const bool inputExtrapolateR,
                            const int inputExtrapolatedRPrecision,
                            const bool inputCorrectNonRadialStrips, const int regionsNumber,
                            const bool doCutOnPrincipals_ = false,
                            const std::string & preEstimatePtDirName = "",
                            const std::string & preEstimateCotThetaDirName = "",
                            const std::string & linearFitLowPtDirName = "",
                            const std::string & linearFitHighPtDirName = "",
                            const std::string & linearFitLongitudinalDirName = "",
                            const bool alignPrincipals = true);

  std::vector<double> estimatedPars() { return estimatedPars_; }
  virtual std::vector<double> principalComponents() { return std::vector<double>(varsNum_, 0.); }
  virtual std::vector<double> normalizedPrincipalComponents() { return std::vector<double>(varsNum_, 0.); }
  double chi2Transverse() const { return chi2Transverse_; }
  int ndofTransverse() const { return ndofTransverse_; }
  double chi2Longitudinal() const { return chi2Longitudinal_; }
  int ndofLongitudinal() const { return ndofLongitudinal_; }
  int ndof() const { return (ndofTransverse_+ndofLongitudinal_); }

 protected:

  virtual bool initialize(const std::vector<double> & vars, const std::vector<int> & layers,
                          const std::vector<int> & stripIndexes) = 0;
  bool consistencyCheck(const std::vector<double> & vars, const std::vector<int> & layers);
  void fillMatrix(std::unordered_map<unsigned long, EstimatorSimple> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  void fillMatrix(std::unordered_map<unsigned long, MatrixReader> * matrices, const unsigned long index,
                  const std::string & fullFileName, const double & deltaPhi, const double & deltaA,
                  const int registerBits, const int bitsPhi, const int bitsA,
                  const int maxBitsMultiplyUnitX, const int maxBitsMultiplyUnitY,
                  const double & scaleFactor, const double & ptSplitValue = 10.,
                  const std::vector<bool> & powerTwoRanges = {});
  std::string buildFullFileName(const std::string & fileName, const std::string & baseDir, const unsigned long & index);

  std::string preEstimatePtDirName_;
  std::string preEstimateCotThetaDirName_;
  std::string linearFitLowPtDirName_;
  std::string linearFitHighPtDirName_;
  std::string linearFitLongitudinalDirName_;
  double ptSplitValue_;
  std::vector<double> principalComponents_;
  std::vector<double> normalizedPrincipalComponents_;
  std::unordered_map<unsigned long, std::vector<double> > meanRadius_;
  unsigned long combinationIndex_;
  std::vector<int> uniqueLayers_;
  unsigned int varsNum_;
  std::string baseDir_;
  CombinationIndexListBuilder combinationIndexListBuilder_;
  bool extrapolateR_;
  int extrapolatedRPrecision_;
  bool correctNonRadialStrips_;
  int regionsNumber_;
  double chi2Transverse_;
  int ndofTransverse_;
  double chi2Longitudinal_;
  int ndofLongitudinal_;
  std::vector<double> estimatedPars_;
  double rotationFactor_;
  bool negativeZ_;
  bool alignPrincipals_;
  bool doCutOnPrincipals_;


  template <class T>
  void fillMatrices(const std::string & baseDir, const std::string & fileName,
                    std::unordered_map<unsigned long, T> * matrices,
                    const double & deltaX = 0., const double & deltaA = 0.,
                    const int registerBits = 0, const int bitsX = 0, const int bitsA = 0,
                    const int maxBitsMultiplyUnitX = 18, const int maxBitsMultiplyUnitY = 27,
                    const double & scaleFactor = 1., const double & ptSplitValue = 10.,
                    const bool normalizeMatrices = false, const bool excludeBarrelFromNormalizedMatrices = false,
                    const std::vector<bool> & powerTwoRanges = {})
  {
    bool fiveOutOfSix = true;

    std::vector<unsigned long> combinationIndexList;
    combinationIndexListBuilder_.fillDefaultIndexList(combinationIndexList, fiveOutOfSix, regionsNumber_);

    for (auto index : combinationIndexList) {
      try {
        std::string fullFileName(fileName);
        fullFileName.replace(fullFileName.find("0"), 1, std::to_string(index));
        fullFileName = baseDir + "/" + fullFileName;
        fillMatrix(matrices, index, fullFileName, deltaX, deltaA, registerBits, bitsX, bitsA,
                   maxBitsMultiplyUnitX, maxBitsMultiplyUnitY, scaleFactor, ptSplitValue,
                   powerTwoRanges);
      }
      catch (int exception) {
        std::cout << "Error: Matrix for combination = " << index << " not found" << std::endl;
        throw;
      }
    }
  }


  /// Align the principal components when there are missing variables
  template <class T>
  void alignPrincipals(std::vector<T> & principals, const int nDof)
  {
    for (int i=0; i<4-nDof; ++i) {
      principals.insert(principals.begin(), 0);
    }
  }
};



#endif //REMOTEPROJECTS_LINEARIZEDTRACKFITTERBASE_H
