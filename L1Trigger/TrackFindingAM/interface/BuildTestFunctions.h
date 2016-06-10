//
// Created by Marco De Mattia on 7/11/15.
//

#ifndef REMOTEPROJECTS_BUILDTESTFUNCTIONS_H
#define REMOTEPROJECTS_BUILDTESTFUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <sstream>
#include "GetVariables.h"


inline void extractCoordinate(const std::vector<double> & vars, const int index, std::vector<double> & coordinate)
{
  for (size_t i=0; i<vars.size()/3; ++i) {
    coordinate.push_back(vars[i*3+index]);
  }
}


void updateMean(std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > & mean,
                const unsigned long combinationIndex, const std::vector<double> & radius);


bool readMean(const std::string & dirName, const std::string & fileName, const unsigned long combinationIndex,
              std::unordered_map<unsigned long, std::vector<double> > & mean);


void initializeVariablesTransformations(const std::vector<std::string> & inputVarNames, const unsigned long combinationIndex,
                                        std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > & variablesTransformations,
                                        const std::string & preEstimateChargeOverPtFileName, const std::string & preEstimateCotThetaFileName,
                                        const std::vector<double> & inputMeanRadius, const std::vector<double> & inputMeanZ);


//template <class T>
//void transformVariables(const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                        const T & variablesTransformations, std::vector<double> & transformedVars,
//                        const double & genChargeOverPt, const double & genCotTheta, const double & genZ0)
//{
//  for (size_t i=0; i<vars.size()/3; ++i) {
//    for (auto v : variablesTransformations) {
//      transformedVars.push_back((*v)(i, vars, uniqueLayers, genChargeOverPt, genCotTheta, genZ0));
//    }
//  }
//}


template <class T>
void transformVariables(const StubsCombination & stubsCombination,
                        const T & variablesTransformations, std::vector<double> & transformedVars)
{
  for (size_t i=0; i<stubsCombination.size(); ++i) {
    for (auto v : variablesTransformations) {
      transformedVars.push_back((*v)(stubsCombination, i));
    }
  }
}


std::vector<std::string> transformedVariablesNames(const int varsSize, const std::vector<std::shared_ptr<TransformBase> > & variablesTransformations);


#endif //REMOTEPROJECTS_BUILDTESTFUNCTIONS_H
