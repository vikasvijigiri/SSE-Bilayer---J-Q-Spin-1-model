import os

def get_matching_directories(root_folder, lamda, j_h, qq_range, l_range):
    matching_directories = []

    for foldername, subfolders, filenames in os.walk(root_folder):
        parts = foldername.split(os.path.sep)

        if len(parts) >= 6:
            lamda_value, j_h_value, qq_value, b_value, l_value = None, None, None, None, None

            for part in parts:
                if part.startswith("lamda_"):
                    lamda_value = part.split('_')[1]
                elif part.startswith("J_H_"):
                    j_h_value = part.split('_')[2]
                    qq_value  = part.split('_')[4]
                    b_value   = part.split('_')[6]
                    #print(j_h_value, qq_value, b_value)
                elif part.startswith("L"):
                    l_value = part[1:]

            #print(qq_value)
            # Check if all values are found and match the specified criteria
            if all(value is not None for value in [lamda_value, j_h_value, qq_value, l_value]):
                if (lamda_value == lamda and j_h_value == j_h and
                        qq_range[0] <= float(qq_value) <= qq_range[1] and
                        l_range[0] <= int(l_value) <= l_range[1]):
                    matching_directories.append(foldername)

    return matching_directories


def submitter(matching_directories):
    # Print the list of matching directories
    for directory_path in matching_directories:
        #print(f"Matching Directory: {directory_path}")
        #print(os.path.join(directory_path, "main"))        
        if os.path.exists(os.path.join(directory_path, "main")):
            os.chdir(directory_path)
            os.system('./job_submitter.sh')  
        else:   
            os.chdir(directory_path)
            os.system('make; make clean')  

# usage
if __name__ == "__main__"
    root_folder = '/home/vikasv/SSE_bilayer/files/Beta_L_by_2'
    lamda_value = '0.0'
    j_h_value = '0.1'
    qq_range = (0.0, 0.17)
    l_range = (4, 16)

    matching_directories = get_matching_directories(root_folder, lamda_value, j_h_value, qq_range, l_range)
    submitter(matching_directories)



