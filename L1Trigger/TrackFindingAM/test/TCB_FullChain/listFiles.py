import os
from os import listdir
from os.path import isfile

def list_files():
    cwd = os.getcwd()
    # The files are in the directory above this one
    files_path = os.path.split(cwd)[0]
    
    only_files = [f for f in listdir(files_path) if f.find("Mu100_PU140_") != -1 and f.find(".root") != -1]
    
    file_list = open("filelist.txt", "w")
    
    for f in only_files:
        # print files_path +"/"+ f
        file_list.write("\"file:"+files_path+"/" + f + "\"\n")
    
    file_list.close()
