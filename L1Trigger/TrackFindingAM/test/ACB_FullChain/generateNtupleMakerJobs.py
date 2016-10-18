import os
from os import listdir
from os.path import isfile

from prepareNtupleMakerBase import prepare_ntuple_maker_cfg

prepare_ntuple_maker_cfg("CB")
prepare_ntuple_maker_cfg("FIT")
prepare_ntuple_maker_cfg("DR")
