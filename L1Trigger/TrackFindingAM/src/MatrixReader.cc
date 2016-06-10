#include "../interface/MatrixReader.h"

using namespace Eigen;

MatrixReader::MatrixReader(const std::string & inputFileName)
{
  // open matrix file and read V and D arrays
  // std::cout << "opening "+inputFileName+" for reading" << std::endl;

  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if (!inputFile) {
    // std::cout << "MatrixReader: Error opening "+inputFileName << std::endl;
    throw 10;
  }

  // Read number of variables and number of track parameters
  inputFile >> nVars_;
  inputFile >> nTrackParameters_;
  nDof_ = nVars_ - nTrackParameters_;
  // std::cout << "Number of variables = " << nVars_ << std::endl;
  // std::cout << "Number of track parameters = " << nTrackParameters_ << std::endl;
  // std::cout << std::endl;

  double x;
  // Read eigenvalues
  sqrtEigenvalues_ = Matrix<long double, Dynamic, 1>::Zero(nVars_);
  for (int i=0; i < nVars_; ++i) {
    inputFile >> x;
    sqrtEigenvalues_(i) = x;
  }
  // std::cout << "sqrt(eigenvalues):" << std::endl;
  // std::cout << sqrtEigenvalues_ << std::endl;

  // Read transformation matrix V from file
  V_ = Matrix<long double, Dynamic, Dynamic>::Zero(nVars_, nVars_);
  for (int i = 0; i < nVars_; ++i) {
    for (int j = 0; j < nVars_; ++j) {
      inputFile >> x;
      V_(i, j) = x;
    }
  }
  // std::cout << "V:" << std::endl;
  // std::cout << std::setprecision(4) << V_ << std::endl;

  meanValues_ = Matrix<long double, Dynamic, 1>::Zero(nVars_);
  for (int i=0; i < nVars_; ++i) {
    inputFile >> x;
    meanValues_(i) = x;
  }
  // std::cout << "meanValues:" << std::endl;
  // std::cout << std::setprecision(4) << meanValues_ << std::endl;

  meanPars_ = Matrix<long double, Dynamic, 1>::Zero(nTrackParameters_);
  for (int i=0; i < nTrackParameters_; ++i) {
    inputFile >> x;
    meanPars_(i) = x;
  }
  // std::cout << "meanTrackParameters:" << std::endl;
  // std::cout << std::setprecision(4) << meanPars_ << std::endl;

  // Read transformation matrix D from file
  D_ = Matrix<long double, Dynamic, Dynamic>::Zero(nTrackParameters_, nVars_);
  for (int i = 0; i < nTrackParameters_; ++i) {
    for (int j = 0; j < nVars_; ++j) {
      inputFile >> x;
      D_(i, j) = x;
    }
  }
  // std::cout << "D:" << std::endl;
  // std::cout << D_ << std::endl;
}


double MatrixReader::normChi2(const Matrix<long double, Dynamic, 1> & vars) const
{
  Matrix<long double, Dynamic, 1> principal = V_*(vars - meanValues_);

  double chi2 = 0.;
  // Use only the constraints to evaluate a chi2
  for (int i=0; i<nDof_; ++i) {
    chi2 += (principal(i)/sqrtEigenvalues_[i])*(principal(i)/sqrtEigenvalues_[i]);
  }
  return chi2/nDof_;
}


std::vector<double> MatrixReader::trackParameters(const Matrix<long double, Dynamic, 1> & vars) const
{
  std::vector<double> pars;

  // Estimate track parameters
  Matrix<long double, Dynamic, 1> estimatedPars = D_ * (vars - meanValues_) + meanPars_;
  for (int i=0; i<nTrackParameters_; ++i) {
    pars.push_back(estimatedPars(i));
  }

  return pars;
}


std::vector<double> MatrixReader::principalComponents(const Matrix<long double, Dynamic, 1> & vars) const
{
  std::vector<double> pcs;

  Matrix<long double, Dynamic, 1> principal = V_*(vars - meanValues_);

  for (int i=0; i<nVars_; ++i) {
    pcs.push_back(principal(i));
  }

  return pcs;
}


std::vector<double> MatrixReader::normalizedPrincipalComponents(const Matrix<long double, Dynamic, 1> & vars) const
{
  std::vector<double> npcs;

  Matrix<long double, Dynamic, 1> principal = V_*(vars - meanValues_);

  for (int i=0; i<nVars_; ++i) {
    npcs.push_back(principal(i)/sqrtEigenvalues_(i));
  }

  return npcs;
}
