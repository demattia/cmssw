import os

index_list = ["2016", "992", "1504", "1760", "1888", "1952", "1984"]

dir_name = "Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10"

for index in index_list:
    os.system("cp ~/d3/Upgrade/FullSimulation/linearizedtrackfit/LinearizedTrackFit/python/ConstantsProduction/"+dir_name+"/*"+index+"*.txt .")
