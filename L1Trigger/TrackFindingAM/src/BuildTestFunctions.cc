//
// Created by Marco De Mattia on 7/11/15.
//

#include "../interface/BuildTestFunctions.h"


void updateMean(std::unordered_map<unsigned long, std::pair<int, std::vector<double> > > & mean,
                const unsigned long combinationIndex, const std::vector<double> & coordinate)
{
  auto it = mean.find(combinationIndex);
  if (it == mean.end()) it = mean.insert(std::make_pair(combinationIndex, std::make_pair(0, std::vector<double>(coordinate.size(), 0)))).first;
  // The counter is common to the whole vector
  it->second.first += 1;
  for (size_t index=0; index<coordinate.size(); ++index) {
    it->second.second.at(index) += (coordinate[index] - it->second.second.at(index)) / it->second.first;
  }
}


bool readMean(const std::string & dirName, const std::string & fileName, const unsigned long combinationIndex,
              std::unordered_map<unsigned long, std::vector<double> > & mean)
{
  auto r = mean.find(combinationIndex);
  if (r == mean.end()) {
    mean.insert(std::make_pair(combinationIndex, std::vector<double>()));
    // Read the radii from the input file
    std::ifstream inputFile;
    inputFile.open(dirName+"/"+fileName+std::to_string(combinationIndex)+".txt");
    if (inputFile) {
      //      std::cout << "reading file: " << dirName+"/MeanRadius_"+std::to_string(combinationIndex)+".txt" << std::endl;
      std::string line;
      std::getline(inputFile, line);
      std::stringstream sline(line);
      double meanValue;
//      std::cout << "mean R = " << std::endl;
      while (sline >> meanValue) {
//        std::cout << meanR << " ";
        mean[combinationIndex].push_back(meanValue);
      }
//      std::cout << std::endl;
    }
    else {
      return false;
    }
  //      else {
  //        std::cout << "readMeanRadius: Error opening "+dirName+"/MeanRadius_"+std::to_string(combinationIndex)+".txt" << std::endl;
  //        // throw;
  //      }
  }
  return true;
}


void initializeVariablesTransformations(const std::vector<std::string> & inputVarNames, const unsigned long combinationIndex,
                                        std::unordered_map<unsigned long, std::vector<std::shared_ptr<TransformBase> > > & variablesTransformations,
                                        const std::string & preEstimateChargeOverPtDirName, const std::string & preEstimateCotThetaDirName,
                                        const std::vector<double> & meanRadius, const std::vector<double> & meanZ)
{
  std::string preEstimateChargeOverPtFileName(preEstimateChargeOverPtDirName+"/matrixVD_"+std::to_string(combinationIndex)+"_pre_chargeOverPt.txt");
  std::string preEstimateCotThetaFileName(preEstimateCotThetaDirName+"/matrixVD_"+std::to_string(combinationIndex)+"_pre_cotTheta.txt");
  std::string preEstimateTgThetaDirName(preEstimateCotThetaDirName);
  int length = preEstimateTgThetaDirName.length();
  // Remove trailing "/" if any
  if (length > 1 && preEstimateTgThetaDirName[length-1] == '/') preEstimateTgThetaDirName = preEstimateTgThetaDirName.substr(0, length - 1);
  std::string preEstimateTgThetaFileName = preEstimateTgThetaDirName+"_tgTheta/matrixVD_"+std::to_string(combinationIndex)+"_pre_tgTheta.txt";

  auto it = variablesTransformations.find(combinationIndex);
  if (it == variablesTransformations.end()) {
    it = variablesTransformations.insert(std::make_pair(combinationIndex, std::vector<std::shared_ptr<TransformBase> >())).first;
    for (auto varName : inputVarNames) {
      if (varName == "phi") {
        variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagatePhi>(varName));
      }
      else if (varName == "CorrectedPhiFirstOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiFirstOrder>(varName, preEstimateChargeOverPtFileName, meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrder>(varName, preEstimateChargeOverPtFileName, meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedR") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedR>(varName, preEstimateChargeOverPtFileName,
                                                                            preEstimateTgThetaFileName, meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder>(varName,
                                                                                       preEstimateChargeOverPtFileName,
                                                                                       preEstimateTgThetaFileName,
                                                                                       meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection>(varName,
                                                                                       preEstimateChargeOverPtFileName,
                                                                                       preEstimateTgThetaFileName,
                                                                                       meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup>(varName,
                                                                                                                     preEstimateChargeOverPtFileName,
                                                                                                                     preEstimateTgThetaFileName,
                                                                                                                     meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_genTheta") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_genTheta>(varName,
                                                                                                                     preEstimateChargeOverPtFileName,
                                                                                                                     preEstimateTgThetaFileName,
                                                                                                                     meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN>(varName,
                                                                                                                         preEstimateChargeOverPtFileName,
                                                                                                                         preEstimateTgThetaFileName,
                                                                                                                         meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN_ExactExtrapolation") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN_ExactExtrapolation>(varName,
                                                                                                                                            preEstimateChargeOverPtFileName,
                                                                                                                                            preEstimateTgThetaFileName,
                                                                                                                                            meanRadius));
      }
      else if (varName == "CorrectedPhiFirstOrderPz") {
        preEstimateChargeOverPtFileName = preEstimateChargeOverPtDirName + "/matrixVD_" + std::to_string(combinationIndex) + "_pre_chargeOverPz.txt";
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiFirstOrderPz>(varName, preEstimateChargeOverPtFileName,
                                                                preEstimateCotThetaFileName, meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderGen") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderGen>(varName, meanRadius));
      }
      else if (varName == "CorrectedPhiSecondOrderGenDeltaZ") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderGenDeltaZ>(varName, meanRadius, meanZ));
      }
      else if (varName == "CorrectedPhiSecondOrderGenExactR") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderGenExactR>(varName, meanRadius, meanZ));
      }
      else if (varName == "CorrectedPhiSecondOrderGenExactRNonRadialStripCorrection") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiSecondOrderGenExactRNonRadialStripCorrection>(varName, meanRadius, meanZ));
      }
      else if (varName == "CorrectedPhiExactGen") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiExactGen>(varName, meanRadius));
      }
      else if (varName == "CorrectedPhiExactGenExactR") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedPhiExactGenExactR>(varName, meanRadius));
      }
      else if (varName == "R") {
        variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagateR>(varName));
      }
      else if (varName == "ExtrapolatedR") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformExtrapolatedR>(varName, preEstimateTgThetaFileName, meanRadius));
      }
      else if (varName == "ExtrapolatedRSecondOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformExtrapolatedRSecondOrder>(varName, preEstimateTgThetaFileName, preEstimateChargeOverPtFileName, meanRadius));
//            std::make_shared<TransformExtrapolatedRSecondOrder>(varName, preEstimateCotThetaFileName, preEstimateChargeOverPtFileName, meanRadius));
      }
      else if (varName == "ExtrapolatedRExact") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformExtrapolatedRExact>(varName));
      }
      else if (varName == "z") {
        variablesTransformations[combinationIndex].push_back(std::make_shared<TransformPropagateZ>(varName));
      }
      else if (varName == "CorrectedZFirstOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedZFirstOrder>(varName, preEstimateCotThetaFileName, meanRadius));
      }
      else if (varName == "CorrectedZSecondOrder") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedZSecondOrder>(varName, preEstimateChargeOverPtFileName,
                                                             preEstimateCotThetaFileName, meanRadius));
      }
      else if (varName == "CorrectedZSecondOrderGen") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedZSecondOrderGen>(varName, meanRadius));
      }
      else if (varName == "CorrectedZExactGen") {
        variablesTransformations[combinationIndex].push_back(
            std::make_shared<TransformCorrectedZExactGen>(varName, meanRadius));
      }
      else {
        std::cout << "Error: variable name " << varName << " not recognized" << std::endl;
        throw;
      }
    }
  }
}


std::vector<std::string> transformedVariablesNames(const int varsSize, const std::vector<std::shared_ptr<TransformBase> > & variablesTransformations)
{
  std::vector<std::string> names;
  for (int i = 0; i < varsSize / 3; ++i) {
    for (auto v : variablesTransformations) {
      names.push_back(v->getName());
    }
  }
  return names;
}
