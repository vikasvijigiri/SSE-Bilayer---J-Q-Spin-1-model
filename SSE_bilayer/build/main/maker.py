import os
import sys

top_dirname = os.getcwd()


Jh = 0.15
fug_values = [0.248, 0.249, 0.25, 0.251, 0.252, 0.253, 
                    0.254,0.255, 0.256, 0.257, 0.258]

#**************** Run this just *************
dirname1 = top_dirname+'/../files/BetaVp_L_by_4/'  
for J_Heisenberg in [Jh]:
    dirname2 = dirname1+'J_H_'+'%.6s' % (str(J_Heisenberg))+'/'
    #for fug in [0.050, 0.080, 0.100, 0.120, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 
    #            0.200, 0.210, 0.215, 0.219, 0.222, 0.223, 0.224, 0.225, 0.226, 0.227, 
    #            0.228, 0.229, 0.230, 0.233, 0.237, 0.242, 0.250, 0.300, 0.400, 0.500]:
    for fug in fug_values:                
                
        dirname3 = dirname2+'fug_'+'%.6s' % (str(fug))+'/'
        for L in [56]:
            dirname4 =  dirname3+'L'+str(L)+'/'
            '''
            num_lines = sum(1 for line in open('output_seeded_data_avg.txt'))
            if num_lines < 15:
            	os.system('sbatch jobscript_dummy.sh')
            '''	
            os.chdir(dirname4)
            os.system('make ')    
            os.system('make clean ')
