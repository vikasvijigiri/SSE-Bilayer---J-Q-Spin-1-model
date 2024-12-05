import os
import sys
#import numpy as np

file_submit = "submit.sh"
node_edit = "node_target.sh"
#os.system("rm -rf " + file_submit)
#kfile = open(file_submit, 'a')


#J_QQ_i = 0.
#J_QQ_f = 1.

#fug = np.array(np.linspace(0.,5.,32))

# np.savetxt('fug.txt', fug, delimiter='\n')   # X is an array

top_dirname = os.getcwd()
for L in [8, 16, 32, 48]:
    dirname1 = top_dirname+'/../files/L'+str(L)+'/'  
    os.system('mkdir -p '+dirname1)
    for J_Heisenberg in [0.0]:
        dirname2 = dirname1+'J_H_'+'%.3s' % (str(J_Heisenberg))+'/'
        os.system('mkdir -p '+dirname2)
        os.system('chmod +x node.sh ')
        for fug in [0.1666]: # close to transition.
            dirname3 = dirname2+'fug_'+'%.5s' % (str(fug))+'/'
            os.system('mkdir -p '+dirname3)            
            os.system('cp node.sh '+dirname3)            
            kfile = open(dirname3 + file_submit, 'w')
            gfile = open(dirname3 + node_edit, 'w')
            #dfile = open(dirname3 + node, 'a')
            #dfile.write('\nsed -g s/sed*/sed -i 14s/.*/#SBATCH --nodelist=bose/ jobscript.sh/g '+node_edit+'\n')
            #dfile.close()             
            #os.system('cp edit_node_target.sh '+dirname3)
            for betaVp in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0]:
                dirname4 =  dirname3+'BetaVp'+str(betaVp)+'/'
                os.system('mkdir -p '+dirname4)
                filename = dirname4+'/input_param.dat'
                pfile = open(filename, 'w')
                #pfile.write('%d %d square-bilayer %f %f 1000000 1000000 10 112' %(L,L,J_Heisenberg,betaVp))
                pfile.write('%d' % (L))
                pfile.write('\n%d' % (1))
                pfile.write('\n%.3f' % (1.))
                pfile.write('\n%.3f' % (fug))
                #pfile.write('\n%.1f' %(J_QQ_f))
                pfile.write('\n%.3f' % (0.))
                pfile.write('\n%.3f' % (betaVp))
                pfile.write('\n10')
                pfile.write('\n100000')
                pfile.write('\n200000')

                pfile.close()
                kfile.write('cd '+dirname4+'\nsbatch jobscript.sh\n')
                gfile.write('cd '+dirname4+'\n'+"sed '14s/.*/#SBATCH --nodelist=bose/g' jobscript_dummy.sh > jobscript.sh"+'\n') 
                #dfile.write('\nsed -g s/#SBATCH --nodelist=bose*/#SBATCH --nodelist=bose/g '+node_edit+'\n')  
                os.system('cp jobscript_dummy.sh '+dirname4)
                os.system('cp Makefile '+dirname4)
                #os.system('cp fug.txt '+dirname4)
            kfile.close()  
            gfile.close()
        '''       
        for J_Biquad in [0.0]:
            dirname3 = dirname2+'J_Bi_'+'%.3s' % (str(J_Biquad))+'/'
            os.system('mkdir -p '+dirname3)
            os.system('cp node.sh '+dirname3)
            kfile = open(dirname3 + file_submit, 'a')
            gfile = open(dirname3 + node_edit, 'a')
            #os.system('cp edit_node_target.sh '+dirname3)                       
            for fug in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]:
                dirname4 = dirname3+'fug_'+'%.3s' % (str(fug))+'/'
                os.system('mkdir -p '+dirname4)
                filename = dirname4+'/input_param.dat'
                pfile = open(filename, 'w')
                #pfile.write('%d %d square-bilayer %f %f 1000000 1000000 10 112' %(L,L,J_Heisenberg,betaVp))
                pfile.write('%d' % (L))
                pfile.write('\n%d' % (1))
                pfile.write('\n%.3f' % (J_Biquad))
                pfile.write('\n%.3f' % (fug))
                #pfile.write('\n%.3f' %(J_QQ_f))
                pfile.write('\n%.3f' % (0.))
                pfile.write('\n%.2f' % (betaVp))
                pfile.write('\n20')
                pfile.write('\n50000')
                pfile.write('\n200000')

                pfile.close()
                kfile.write('cd '+dirname4+'\nsbatch jobscript.sh\n')
                gfile.write('cd '+dirname4+'\n'+"sed -i '14s/.*/#SBATCH --nodelist=bose1/g' jobscript.sh"+'\n')
                os.system('cp jobscript.sh '+dirname4)
                os.system('cp Makefile '+dirname4)
            kfile.close()
            gfile.close()
            #os.system('cp fug.txt '+dirname4)
        '''    
                

#kfile.close()
