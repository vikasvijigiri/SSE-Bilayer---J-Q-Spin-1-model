import os
import pandas as pd
import re
import numpy as np



# **********************************************************************************************************
#              GENERIC PART (NEVER CHANGE (except what you need exactly need is not found) 
# **********************************************************************************************************


class depth_maker:
    def __init__(self, mdir):
        # System parameters
        self.mdir = mdir
        print("Hello, you are extracting *.txt files from ", mdir)
        


    def extract_dat_files(self, root_folder, endswith):
        # Create the output folder if it doesn't exist
        # Create the output folder if it doesn't exist
        #os.makedirs(output_folder, exist_ok=True)

        #os.chmod(output_folder, 0o777)

        found_files_list = []
        # Iterate through all subdirectories and extract *.dat files
        for foldername, subfolders, filenames in os.walk(root_folder):
            #print(foldername,subfolders,filenames)
            for filename in filenames:
                if filename.endswith(endswith):
                    file_path = os.path.join(foldername, filename)
                    #print(f"Found: {file_path}")
                    depth = foldername.count(os.path.sep) - root_folder.count(os.path.sep)

                    found_files_list.append(file_path)                    

                    # Copy the *.dat file to the output folder
                    #shutil.copy(file_path, output_folder)
                    #print(f"Extracted: {file_path} to {output_folder}/{filename}")
        return found_files_list

    def print_found_files(self, found_files_list): 
        # Now found_files_list contains all the file paths with the specified extension
        print("List of found files:")
        for file_path in found_files_list:
            print(file_path)


# *********************************************************************************************
#                             MAIN (Access)
# *********************************************************************************************



def run_depth_maker(root_folder, endswith="data_avg.txt"):  
    dmak = depth_maker(root_folder)      
    found_files_list = dmak.extract_dat_files(root_folder, endswith)
    dmak.print_found_files(found_files_list)
    for file_path in found_files_list:
        directory_path = os.path.dirname(file_path)
        if not os.path.exists(os.path.join(directory_path, "main")):
            os.chdir(directory_path)
            os.system('make ')    
            os.system('make clean ')



if __name__ == "__main__":
    root_folder = "/home/vikasv/SSE_bilayer/files/Beta_L_by_2/"   # change this to your folder you want.
    run_depth_maker(root_folder, "input_param.dat")
