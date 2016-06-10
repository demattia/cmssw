#ifndef MATRIXREADER_H
#define MATRIXREADER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "Eigen/Eigenvalues"
//#include <Eigenvalues>

using namespace Eigen;

class MatrixReader
{
public:
  MatrixReader(const std::string & inputFileName);
  double normChi2(const Matrix<long double, Dynamic, 1> & vars) const;
  std::vector<double> trackParameters(const Matrix<long double, Dynamic, 1> & vars) const;
  std::vector<double> principalComponents(const Matrix<long double, Dynamic, 1> & vars) const;
  std::vector<double> normalizedPrincipalComponents(const Matrix<long double, Dynamic, 1> & vars) const;
  int nDof() { return nDof_; }

private:
  int nVars_;
  int nTrackParameters_;
  int nDof_;
  Matrix<long double, Dynamic, 1> sqrtEigenvalues_;
  Matrix<long double, Dynamic, 1> meanValues_;
  Matrix<long double, Dynamic, 1> meanPars_;
  Matrix<long double, Dynamic, Dynamic> V_;
  Matrix<long double, Dynamic, Dynamic> D_;
};

#endif // MATRIXREADER_H
