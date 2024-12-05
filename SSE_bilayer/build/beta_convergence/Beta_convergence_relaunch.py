import os
import sys
import numpy as np



top_dirname = os.getcwd() + "/../files/Beta_convergence/"
for J_Heisenberg in [0.025]:
    dirname1 = top_dirname + 'J_H_'+'%.5s' % (str(J_Heisenberg))+'/'
    os.system('mkdir -p '+dirname1) 
    for lmbda in [0, 0.5, 1.0]:
        dirname2 = dirname1 + 'lamda_'+'%.5s' % (str(lmbda))+'/'
        os.system('mkdir -p '+dirname2)
        for fug in [0.160,0.180,0.210]: 			 			
            dirname3 = dirname2 + 'fug_'+'%.5s' % (str(fug))+'/'
            os.system('mkdir -p '+dirname3)            
            for beta in [0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 4.0, 6.0, 8.0]:
                dirname4 = dirname3 + 'Beta_'+str(beta)+'/'  
                os.system('mkdir -p '+dirname1)
                for L in [4]: #8,16,32]:
                    dirname5 =  dirname4+'L'+str(L)+'/'
                    os.system('mkdir -p '+dirname5)

                    filename = dirname5+'/input_param.dat'
                    pfile = open(filename, 'w')
                    pfile.write('%d' % (L))
                    pfile.write('\n%d' % (2))
                    pfile.write('\n%.3f' % (J_Heisenberg))
                    pfile.write('\n%.3f' % (fug))
                    pfile.write('\n%.3f' % (1.))
                    pfile.write('\n%.3f' % (beta))
                    pfile.write('\n10')
                    pfile.write('\n200000')
                    pfile.write('\n400000')
                    pfile.write('\n1.')
                    pfile.write('\n0')

                    pfile.close()
                    os.system('cp jobscript_dummy.sh '+dirname5)
                    os.system('cp Makefile '+dirname5)
