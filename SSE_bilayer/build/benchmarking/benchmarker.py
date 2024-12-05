import os
import sys
import numpy as np


JQBvals = [[1.0,1.0,0.2], [0.9,0.4,0.3], [0.6,0.5,0.5], [0.3,0.8,0.2], [1.15,0.88,0.12]]  
lambda_vals = [1., 0.5, 0., 0.73, 0.4]

top_dirname = os.getcwd() + "/../../files/benchmarking/"
for ii, JQB in enumerate(JQBvals):
    JH, JQQ, JB = JQB 
    dirname1 = top_dirname + 'J_H_'+'%.5s' % (str(JH)) + '_QQ_'+'%.5s' % (str(JQQ)) + '_B_'+'%.5s' % (str(JB)) +'/'
    os.system('mkdir -p '+dirname1) 
    lmbda = lambda_vals[ii]
    dirname2 = dirname1 + 'lamda_'+'%.5s' % (str(lmbda))+'/'
    os.system('mkdir -p '+dirname2)
    for beta in [20]:
        dirname3 = dirname2 + 'Beta_'+str(beta)+'/'  
        os.system('mkdir -p '+dirname1)
        for L in [4]: #8,16,32]:
            dirname4 =  dirname3+'L'+str(L)+'/'
            os.system('mkdir -p '+dirname4)

            filename = dirname4+'/input_param.dat'
            pfile = open(filename, 'w')
            pfile.write('%d' % (L))
            pfile.write('\n%d' % (2))
            pfile.write('\n%.3f' % (JH))
            pfile.write('\n%.3f' % (JQQ))
            pfile.write('\n%.3f' % (JB))
            pfile.write('\n%.3f' % (beta))
            pfile.write('\n200000')
            pfile.write('\n400000')
            pfile.write('\n%.5f' % (lmbda))
            pfile.write('\n0')

            pfile.close()
            os.system('cp jobscript_dummy.sh '+dirname4)
            os.system('cp Makefile '+dirname4)
            os.system('cp job_submitter.sh '+dirname4)
