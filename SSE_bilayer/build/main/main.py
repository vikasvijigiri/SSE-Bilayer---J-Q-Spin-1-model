import os
import sys
import numpy as np


# Create grids and submit jobs

def generate_custom_distribution(num_points, mean, left_boundary, right_boundary):
    # Calculate the spacing factors based on the boundaries
    left_spacing_factor = 1. * (mean - left_boundary)
    right_spacing_factor = 1. * (right_boundary - mean)

    # Generate 1D grid data with quadratic spacing
    left_spacing = -left_spacing_factor * (np.linspace(0, 1, num_points // 2))**1.8
    right_spacing = right_spacing_factor * (np.linspace(0, 1, num_points // 2))**1.7

    # Combine left and right spacing
    spacing = np.concatenate((left_spacing[::-1], right_spacing[1:]))
    grid_data = mean + spacing

    # Clip the values to ensure they are within the specified range [left_boundary, right_boundary]
    grid_data = np.clip(grid_data, left_boundary, right_boundary)
    #print(grid_data)
    # Round off grid_data to three decimals
    grid_data = np.round(grid_data, 3)

    return grid_data


# Set the parameters
JH = 0.1
JB = 1.
lambda_vals = [0.0]


num_points = 24
gc = 0.083
left_boundary = 0.0
right_boundary = 0.17 
# Generate custom distribution
fug = generate_custom_distribution(num_points, gc, left_boundary, right_boundary)
#[0., 0.01, 0.25, 0.04, 0.05, 0.07, 0.1]


top_dirname = os.getcwd() + "/../../files/Beta_L_by_2/"
for lmbda in lambda_vals:
    dirname1 = top_dirname + 'lamda_'+'%.5s' % (str(lmbda))+'/'
    os.system('mkdir -p '+dirname1)
    for ii, JQQ in enumerate(fug):
        dirname2 = dirname1 + 'J_H_'+'%.5s' % (str(JH)) + '_QQ_'+'%.5s' % (str(JQQ)) + '_B_'+'%.5s' % (str(JB)) +'/'
        os.system('mkdir -p '+dirname2) 
        for L in [4, 8, 12, 16, 20, 24, 32, 40]: #8,16,32]:
            dirname3 =  dirname2+'L'+str(L)+'/'
            os.system('mkdir -p '+dirname3)

            filename = dirname3+'/input_param.dat'
            pfile = open(filename, 'w')
            pfile.write('%d' % (L))
            pfile.write('\n%d' % (L))
            pfile.write('\n%.3f' % (JH))
            pfile.write('\n%.3f' % (JQQ))
            pfile.write('\n%.3f' % (JB))
            pfile.write('\n%.3f' % (L/2))
            pfile.write('\n200000')
            pfile.write('\n400000')
            pfile.write('\n%.5f' % (lmbda))
            pfile.write('\n0')

            pfile.close()
            os.system('cp jobscript_dummy.sh '+dirname3)
            os.system('cp Makefile '+dirname3)
            os.system('cp job_submitter.sh '+dirname3)
