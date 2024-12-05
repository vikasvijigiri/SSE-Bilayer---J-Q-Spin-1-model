#include <iostream>

#include "../include/SSEvariables.hpp"
#include "../include/SSElattice.hpp"
#include "../include/SSEobservables.hpp"
#include "../include/SSEupdates.hpp"



static ran rann;
using namespace std;


double computeWeight(SSEvariables& , SSElattice* , double , int , int& ); 
void applyOperator(SSEvariables& , SSElattice* , int* );
void applyOperator_measurement(SSEvariables& , SSElattice* , int* , int& , int& );

// ****************************************************************************************
//                                  Main MC parts 
// ****************************************************************************************

void
SSEupdates::mcstep(SSEvariables& sv, SSElattice* lattice) {
  diag_update(sv, lattice);
  looper(sv, lattice);
  checkl(sv);
}

void
SSEupdates::mcstep_measurement(SSEvariables& sv, SSElattice* lattice, SSEobservables& obs) {
  diag_update_measurement(sv, lattice, obs);
  looper_measurement(sv, lattice, obs);
}


// ****************************************************************************************
//                                  Updates part
// ****************************************************************************************


// Local helpful function
double computeWeight(SSEvariables& sv, SSElattice* lattice, double awgt[][2], int k, int& b) {
    double p = 1.;
    int kk = static_cast <int> (k/4);
    if (kk == 0) {
        b = (k%4) + 4*int(rann() * 2*sv.Ns);
        int ss1 = (1 + lattice[sv.JHsites[b][0]].S()) / 2;
        int ss2 = (1 + lattice[sv.JHsites[b][1]].S()) / 2;
        p *= awgt[ss1][ss2];
        //p *= sv.Hsgn[b];
    } else if (kk == 1) {
        b = (k%4) + 4*int(rann() * 2*sv.Ns);
        //b = k*sv.Ns + int(rann() * sv.NQ); // Old method
        for (int j = 0; j < 4; ++j) {
            int ss1 = (1 + lattice[sv.JQsites[b][2 * j]].S()) / 2;
            int ss2 = (1 + lattice[sv.JQsites[b][2 * j + 1]].S()) / 2;
            p *= awgt[ss1][ss2];
        }
        //p *= sv.QQsgn[b];
    } else {
        b = (k%4) + 2*int(rann() * 2*sv.Ns);        
        //b = k*sv.Ns + int(rann() * sv.NB);
        for (int j = 0; j < 2; ++j) {
            int ss1 = (1 + lattice[sv.JBsites[b][2 * j]].S()) / 2;
            int ss2 = (1 + lattice[sv.JBsites[b][2 * j + 1]].S()) / 2;
            p *= awgt[ss1][ss2];
        }
        //p *= sv.Bsgn[b];
    }
    return p;
}


// Local helpful function
void applyOperator(SSEvariables& sv, SSElattice* lattice, int** ostr, int i) {
    int b = ostr[3][i] - 1;
    int t = ostr[2][i];

    
    if (t == 0) {
        int op = ostr[0][i];
        if (op == 2) {
            lattice[sv.JHsites[b][0]].flip();
            lattice[sv.JHsites[b][1]].flip();
        }    
    } else if (t == 1) {
        for (int j = 0; j < 2; j++) {
            int op = ostr[j][i];
            int index = 4 * j;
//            if (op > 4 || op < 0) {
//                //std::cout << op << std::endl;
//                printf("wrong operator QQ!! ");
//                exit(0);
//            }
            if (op == 2 || op == 4) {
                lattice[sv.JQsites[b][0 + index]].flip();
                lattice[sv.JQsites[b][1 + index]].flip();
            }
            if (op == 3 || op == 4) {
                lattice[sv.JQsites[b][2 + index]].flip();
                lattice[sv.JQsites[b][3 + index]].flip();
            }
        }
    } else if (t == 2) {
        int op = ostr[0][i];
//        if (op > 4 || op < 0) {
//            //std::cout << op << std::endl;
//            printf("wrong operator BB!! ");
//            exit(0);
//        }
        if (op == 2 || op == 4) {
            lattice[sv.JBsites[b][0]].flip();
            lattice[sv.JBsites[b][1]].flip();
        }
        if (op == 3 || op == 4) {
            lattice[sv.JBsites[b][2]].flip();
            lattice[sv.JBsites[b][3]].flip();
        }
    }
//    else {
//        //std::cout << t << std::endl;
//        printf("wrong operator!! ");
//        exit(0);
//    }
}



// Diagonal updates
void
SSEupdates::diag_update(SSEvariables& sv, SSElattice* lattice) {
  double p, r, cp;
  int b;

//  int dum_latt[32]; // for a 4 x 4 lattice

//  for (int i = 0; i < sv.Nm; ++i) {
//      dum_latt[i] = lattice[i].S();
//  }

  for (long i = 0; i < sv.Lc; ++i) {
    int o = str[0][i] * str[1][i];
    if (o == 0) {     // Encountered identity operator, try to insert a diagonal operator
      r = rann();
      cp = 0.;
      for (int kk = 0; kk < 10; kk++) {
        int k = static_cast <int> (kk/4);
        cp += sv.cum_prob[kk];
        if (cp > r && sv.cum_prob[kk] > 1e-6) { // J or Q or Biquad?
          p = computeWeight(sv, lattice, awgt, kk, b);
          p = p * sv.prob_in / static_cast<double>(sv.Lc - sv.n1);
          if (rann() <= p) 
          {
            str[0][i] = 1;      // Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to QoQ and Bi and H);  1 for diag 2 for off-diag
            str[1][i] = 1;      // Off-digonal (or On-diag) operator within the 2nd biquad bond (relevant to QoQ));  1 for diag 2 for off-diag
            str[2][i] = k;      // Which type of operator, 0 for H, 1 for  QQ, 2 for Biquad;
            str[3][i] = b + 1;  // Bond-index;
            ++Nd[k];
            ++sv.n1;
            //break;
          }
          break;
        }
      }
    } else if (o == 1) {    // Encountered diagonal operators, attempt to remove
      p = sv.prob_rm * static_cast<double>(sv.Lc - sv.n1 + 1);
      if (rann() <= p) {
        --Nd[str[2][i]];
        --sv.n1;
        for (int j = 0; j < 4; ++j) str[j][i] = 0;
      }
    } 
    else {  // Encountered offdiagonal operators
        applyOperator(sv, lattice, str, i); 
    }
  }


//  for (int i = 0; i < sv.Nm; ++i) {
//      if (dum_latt[i] != lattice[i].S()) {
//        printf("Some problem");
//        exit(0);         
//      }
//  }
		
}


// ********************************** Equilibration Looper update *****************************************

// offdiagonal/looper updates
void
SSEupdates::looper(SSEvariables& sv, SSElattice * lattice) {

  long a = (Nd[0] + 4 * Nd[1] + 2 * Nd[2]);
  long * X = new long[4 * a];
  long * lpos = new  long[a];
  int * visited_vertex = new int[sv.Lc];
  
  // Linked list construction
  long v0 = 0;
  for (long i = 0; i < sv.Lc; ++i) {
      int o = str[0][i];
      visited_vertex[i] = 0;
      if (o != 0) {
          int b = str[3][i] - 1;
          int v = str[2][i];
          for (int j = 0; j < type[v]; j++) {
              int s1, s2;
              if (v == 0) {
                s1 = sv.JHsites[b][0];
                s2 = sv.JHsites[b][1];
              } else if (v == 1) {
                s1 = sv.JQsites[b][2 * j];
                s2 = sv.JQsites[b][2 * j + 1];
              } else if (v == 2) {
                s1 = sv.JBsites[b][2 * j];
                s2 = sv.JBsites[b][2 * j + 1];
              } else {
                printf("wrong vertex ");
                exit(0);
              }
              lpos[static_cast<long>(v0 / 4)] = 4 * i + j;
              long v1 = last[s1];
              long v2 = last[s2];

              if (v1 != -1) {
                  X[v1] = v0;
                  X[v0] = v1;
              } else {
                  frst[s1] = v0;
              }
              if (v2 != -1) {
                  X[v2] = v0 + 1;
                  X[v0 + 1] = v2;
              } else {
                  frst[s2] = v0 + 1;
              }
              last[s1] = v0 + 2;
              last[s2] = v0 + 3;
              v0 += 4;
          }
      }
  }

  
  // PBC loops 
  for (int k = 0; k < sv.Nm; ++k) {
    long v1 = frst[k];
    if (v1 != -1) {
      long v2 = last[k];
      X[v2] = v1;
      X[v1] = v2;
    }
  }


  sv.nl = 10*sv.lx;
  // Main looper 
  if (v0 != 0) {
      for (int j = 0; j < sv.nl; ++j) {
          long vi = static_cast<long>(rann() * v0);
          long v1 = vi;
          if (X[vi] < 0) continue;
          while (true) {
              long l = lpos[static_cast<long>(v1 / 4)];
              long i = static_cast<long>(l / 4);
              int k = str[2][i];
              if (k == 0) {
                  str[0][i] = 3 - str[0][i];
              } 
              else if (k == 1) {
                  int v = static_cast<int>((l % 4) / 2);        // Index of which bond out of those 4
                  str[v][i] = flip[str[v][i]][(l % 4) % 2];
              } 
              else if (k == 2) {
                  int v = static_cast<int>((l % 4) % 2);
                  str[0][i] = flip[str[0][i]][v];
              }
              long v2 = v1 ^ 1;
              X[v1] = -1;
              v1 = X[v2];
              X[v2] = -1;
              visited_vertex[i] += 1;
              if (v1 == vi || v2 == vi) break;
          }
      }
  }

//  // Update number of loops, nl
//  int number_visited = 0;
//  for (int j = 0; j < sv.Lc; j++) {
//    //i = int(lpos[j] / Np);
//    int k = str[2][j];
//    number_visited += visited_vertex[j] / type[k];
//  }
//	
//  if (number_visited  < sv.n1) {
//    sv.nl += int(sv.nl / 3);
//  }
//  
 
  // Update spin configuration here
  for (int j = 0; j < sv.Nm; ++j) {
    if (last[j] != -1) {
      if (X[last[j]] == -1){
          lattice[j].flip();
      }
    } else {
      if (rann() < 0.5) lattice[j].flip();
    }
    last[j] = -1;
    frst[j] = -1;
  }

  delete[] lpos;
  delete[] X;
  delete[] visited_vertex;
}


// ****************************************************************************************
//                UPdates with measurements
// ****************************************************************************************


void applyOperator_measurement(SSEvariables& sv, SSElattice* lattice, int** ostr, int i, int& jjx, int& jjy) {
    int b = ostr[3][i] - 1;
    int t = ostr[2][i];

    if (t == 0) {
        int op = ostr[0][i];
        if (op == 2){
            lattice[sv.JHsites[b][0]].flip();
            lattice[sv.JHsites[b][1]].flip();
				    jjx += sv.JHsgnx[b]*lattice[sv.JHsites[b][1]].S();
				    jjy += sv.JHsgny[b]*lattice[sv.JHsites[b][1]].S();
        }
    } else if (t == 1) {
        for (int j = 0; j < 2; j++) {
            int op = ostr[j][i];
            int index = 4 * j;
            if (op == 2 || op == 4) {
                lattice[sv.JQsites[b][0 + index]].flip();
                lattice[sv.JQsites[b][1 + index]].flip();
				        jjx += sv.JQsgnx[b]*lattice[sv.JQsites[b][1 + index]].S();
				        jjy += sv.JQsgny[b]*lattice[sv.JQsites[b][1 + index]].S();
            }
            if (op == 3 || op == 4) {
                lattice[sv.JQsites[b][2 + index]].flip();
                lattice[sv.JQsites[b][3 + index]].flip();
				        jjx += sv.JQsgnx[b]*lattice[sv.JQsites[b][3 + index]].S();
				        jjy += sv.JQsgny[b]*lattice[sv.JQsites[b][3 + index]].S();
            }
        }
    } else if (t == 2) {
        int op = ostr[0][i];
        if (op == 2 || op == 4) {
            lattice[sv.JBsites[b][0]].flip();
            lattice[sv.JBsites[b][1]].flip();
				    jjx += sv.JBsgnx[b]*lattice[sv.JBsites[b][1]].S();
            jjy += sv.JBsgny[b]*lattice[sv.JBsites[b][1]].S();
        }
        if (op == 3 || op == 4) {
            lattice[sv.JBsites[b][2]].flip();
            lattice[sv.JBsites[b][3]].flip();
				    jjx += sv.JBsgnx[b]*lattice[sv.JBsites[b][3]].S();
            jjy += sv.JBsgny[b]*lattice[sv.JBsites[b][3]].S();
        }
    }
}


// Diagonal update with measurements
void
SSEupdates::diag_update_measurement(SSEvariables& sv, SSElattice* lattice, SSEobservables& obs) {
  double p, r, cp;
  int b;
  int jjx=0, jjy=0;

  // Measurement of some diagonal quantities.
  double mzs1 = 0., mzs2 = 0.;  
  for(int x = 0; x < sv.Ns; x++){
      double f1 = 0.25 * lattice[0].S() * lattice[x].S();                             // Layer 1
      double f2 = 0.25 * lattice[sv.Ns].S() * lattice[sv.Ns+x].S();                   // Layer 2
      double f3 = 0.25 * lattice[sv.Ns].S() * lattice[x].S();                         // inter layers
      double f4 = 0.25 * lattice[0].S() * lattice[sv.Ns+x].S();                       // inter layers
	
      double g1 = 0.5 * lattice[x].S()  * pow(-1, (x / sv.lx) + (x % sv.lx));         // Layer 1
      double g2 = 0.5 * lattice[x+sv.Ns].S() * pow(-1, (x / sv.lx) + (x % sv.lx));		// Layer 2
      obs.B[0][x] += f1;  // Layer 1
			obs.B[1][x] += f2;  // Layer 2
			obs.B[2][x] += f3;  // Inter layers
      obs.B[3][x] += f4;  // Inter layers       
			mzs1 += g1;
			mzs2 += g2;	
	}

	mzs1 /= static_cast<double>(sv.Nm);
  mzs2 /= static_cast<double>(sv.Nm);

  // ~ <Magnetization>
  obs.SMagL1  += mzs1;                          // Layer 1
  obs.SMagL2  += mzs2;                          // Layer 2
  
  // ~ <Magnetization^2>
  obs.SMag_squareL1 += mzs1*mzs1;               // Layer 1
  obs.SMag_squareL2 += mzs2*mzs2;               // Layer 2
  obs.SMag_squareI1 += mzs1*mzs2 + mzs2*mzs1;    // Inter layer

  // ~ <Magnetization^4>
  obs.SMag_fourL1 += pow(mzs1,4);                     // Layer 1 
  obs.SMag_fourL2 += pow(mzs2,4);                     // Layer 2
  obs.SMag_fourI1 += 6.*pow(mzs1,2)*pow(mzs2,2);      // Inter layer
  obs.SMag_fourI2 += 4.*mzs1*pow(mzs2, 3);            // Inter layer
  obs.SMag_fourI3 += 4.*mzs2*pow(mzs1, 3);            // Inter layer 

  // update with stiffness measurements;
  for (long i = 0; i < sv.Lc; ++i) {
    int o = str[0][i] * str[1][i];
    if (o == 0) {
      r = rann();
      cp = 0.;
      for (int kk = 0; kk < 10; kk++) {
        int k = static_cast <int> (kk/4);
        cp += sv.cum_prob[kk];
        if (cp > r && sv.cum_prob[kk] > 1e-6) { // J or Q or Biquad?
          p = computeWeight(sv, lattice, awgt, kk, b);
          p = p * sv.prob_in / static_cast<double>(sv.Lc - sv.n1);
          if (rann() <= p) 
          {
            str[0][i] = 1;      // Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to QoQ and Bi and H);  1 for diag 2 for off-diag
            str[1][i] = 1;      // Off-digonal (or On-diag) operator within the 2nd biquad bond (relevant to QoQ));  1 for diag 2 for off-diag
            str[2][i] = k;      // Which type of operator, 0 for H, 1 for  QQ, 2 for Biquad;
            str[3][i] = b + 1;  // Bond-index;
            ++Nd[k];
            ++sv.n1;
            //break;          
          }
         break;
        }
      }
    } else if (o == 1) {
      p = sv.prob_rm * static_cast<double>(sv.Lc - sv.n1 + 1);
      if (rann() <= p) {
        --Nd[str[2][i]];
        --sv.n1;
        for (int j = 0; j < 4; ++j) str[j][i] = 0;
      }
    } else {
        applyOperator_measurement(sv, lattice, str, i, jjx, jjy);
    }
  }	
 

  // Avg Stiffness along x and y;
  obs.stiffx += 0.5 * static_cast<double>(pow(jjx, 2)) / static_cast<double>(sv.lx*sv.lx)/sv.Beta;    // Layer 1 & 2; Note: you cant seperate them out.
  obs.stiffy += 0.5 * static_cast<double>(pow(jjy, 2)) / static_cast<double>(sv.ly*sv.ly)/sv.Beta;    // Layer 1 & 2; Note: you cant seperate them out.


  double e = -static_cast<double>(sv.n1)/static_cast<double>(sv.Nm)/sv.Beta;

  // ~ <Energy> 
  obs.enrg += e;	

  // ~ <Energy^2> 
  obs.enrg2 += pow(e, 2);	

  // ~ <Energy^4> 
  obs.enrg4 += pow(e, 4);	
}




// offdiagonal/looper update with measurements
void
SSEupdates::looper_measurement(SSEvariables& sv, SSElattice* lattice, SSEobservables& obs) {

  long a = (Nd[0] + 4 * Nd[1] + 2 * Nd[2]);
  long* X = new long[4 * a];
  long* lpos = new  long[a];
  bool* flag = new bool [sv.Lc];


  
  // Linked list construction
  long v0 = 0;
  for (long i = 0; i < sv.Lc; ++i) {
      int o = str[0][i]; flag[i] = 0;
      if (o != 0) {
          int b = str[3][i] - 1;
          int v = str[2][i];
          for (int j = 0; j < type[v]; j++) {
              int s1, s2;
              if (v == 0) {
                s1 = sv.JHsites[b][0];
                s2 = sv.JHsites[b][1];
              } else if (v == 1) {
                s1 = sv.JQsites[b][2 * j];
                s2 = sv.JQsites[b][2 * j + 1];
              } else if (v == 2) {
                s1 = sv.JBsites[b][2 * j];
                s2 = sv.JBsites[b][2 * j + 1];
              } else {
                printf("wrong vertex ");
                exit(0);
              }
              lpos[static_cast<long>(v0 / 4)] = 4 * i + j;
              long v1 = last[s1];
              long v2 = last[s2];

              if (v1 != -1) {
                  X[v1] = v0;
                  X[v0] = v1;
              } else {
                  frst[s1] = v0;
              }
              if (v2 != -1) {
                  X[v2] = v0 + 1;
                  X[v0 + 1] = v2;
              } else {
                  frst[s2] = v0 + 1;
              }
              last[s1] = v0 + 2;
              last[s2] = v0 + 3;
              v0 += 4;
          }
      }
  }

  


  // PBC loops 
  for (int k = 0; k < sv.Nm; ++k) {
    long v1 = frst[k];
    if (v1 != -1) {
      long v2 = last[k];
      X[v2] = v1;
      X[v1] = v2;
    }
  }




  // ****************** Measurements ****************** 



  // Correlators here
  for (long jj = 0; jj < v0/4; ++jj) {
  		long ii = jj*4;
      long l = lpos[jj];
      long i = long(l / 4);
      if (flag[i] == 0) {
        int o = str[0][i] * str[1][i];
        flag[i] = 1;
        if (o != 0) {
          int k = str[2][i];
          //<H-Bond.H-Bond>
          if (k == 0) { // k = 0 --> H,   1--> QoQ,   2--> Bi
            int b = str[3][i] - 1;
            int vb2 = int(b / (4*sv.Ns));
            // Bond indices
            int b3 = b; //long(b / 4);                     // x-bond
            int b5 = b % (4*sv.Ns);//int((b % (4*sv.Ns)) / 4);       // y-bond
//            int b3l = b % 4;
//            int b5l = (b % (4*sv.Ns)) % 4;
  
            l = lpos[long(((ii + 4) % v0) / 4)];
            int j = long(l / 4) % sv.Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 0) { // k = 0 --> H,   1--> QoQ,   2--> Bi
                b = str[3][j] - 1;
                int vb1 = long(b / ( 4*sv.Ns ));
                // Bond indices
                int b2 = b; //long(b / 4);                 // x-bond
                int b4 = b % (4*sv.Ns);  //long((b % ( 4*sv.Ns )) / 4); // y-bond
//                int b2l = b % 4;
//                int b4l = (b % (4*sv.Ns)) % 4;
                if (vb1 == 0 and vb2 == 0) {  // xx-type; So, b2, b3 are horizontal bonds
                  obs.Cx[b2][b3] += sv.n1-1;
                  //if (b2 != b3) obs.Cx[b3][b2] += sv.n1-1;
                }
                if (vb1 == 1 and vb2 == 1) {  // yy-type; So, b4, b5 are vertical bonds
                  obs.Cy[b5][b4] += sv.n1-1;
                  //if (b5 != b4) obs.Cy[b4][b5] += sv.n1-1;
                }
              }
            }
          }
          //<Biquad.Biquad>
          else if (k == 2) { // k = 0 --> H,   1--> QoQ,    2--> Bi
            int b = str[3][i] - 1;
            int vb2 = long(b / (2*sv.Ns ));
            // Bond indices
            int b3 = b;//long(b / 2);                       // x-bond
            int b5 = (b % (2*sv.Ns ));//long((b % (2*sv.Ns )) / 2);        // y-bond
            //int b3l = b % 2;
            //int b5l = (b % (2*sv.Ns )) % 2;

            l = lpos[long(((ii + 8) % v0) / 4)];
            int j = long(l / 4) % sv.Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 2) { // k = 0 --> H,   1--> QoQ,    2--> Bi
                b = str[3][j] - 1;
                int vb1 = long(b / (2*sv.Ns ));
                // Bond indices
                int b2 = b;//long(b / 2);                   // x-bond
                int b4 = (b % (2*sv.Ns )); // long((b % (2*sv.Ns )) / 2);    // y-bond
                //int b2l = b % 2;
                //int b4l = (b % (2*sv.Ns )) % 2;

                if (vb1 == 0 and vb2 == 0) { 	// xx type. So, b2, b3 are horizontal bonds
                  obs.Dx[b2][b3] += sv.n1-1;
                  //obs.Dx[b3][b2] += sv.n1-1;
                  //if (b2 != b3) obs.Dx[b3][b2] += sv.n1-1;

//                  obs.Ex[b2/2][b3/2] += sv.n1-1;
//                  //obs.Dx[b3][b2] += sv.n1-1;
//                  if (b2 != b3) obs.Ex[b3/2][b2/2] += sv.n1-1;

                }
                if (vb1 == 1 and vb2 == 1) {	// yy type. So, b4, b5 are vertical bonds
                  obs.Dy[b4][b5] += sv.n1-1;
                  //obs.Dy[b5][b4] += sv.n1-1;
                  //if (b4 != b5) obs.Dy[b5][b4] += sv.n1-1;

//                  obs.Ey[b4/2][b5/2] += sv.n1-1;
//                  //obs.Dy[b5][b4] += sv.n1-1;
//                  if (b4 != b5) obs.Ey[b5/2][b4/2] += sv.n1-1;
                }
              }
            }
          }
        }
      }
  }

  // Back to updates part
  // Main looper 
  for (int j = 0; j < sv.nl; ++j) {
      long vi = static_cast<long>(rann() * v0);
      long v1 = vi;
      if (X[vi] < 0) continue;
      while (true) {
          long l = lpos[static_cast<long>(v1 / 4)];
          long i = static_cast<long>(l / 4);
          //int ic = static_cast<int>(v1 % 4);
          int k = str[2][i];
          if (k == 0) {
              str[0][i] = 3 - str[0][i];
          } else if (k == 1) {
              int v = static_cast<int>((l % 4) / 2); // Index of which bond out of those 4
              str[v][i] = flip[str[v][i]][(l % 4) % 2];
          } else if (k == 2) {
              int v = static_cast<int>((l % 4) % 2);
              str[0][i] = flip[str[0][i]][v];
          }
          //long v2 = 4 * (static_cast<long>(v1 / 4)) + oc;
          X[v1] = -1;
          long v2 = v1 ^ 1;
          v1 = X[v2];
          X[v2] = -1;
          if (v1 == vi || v2 == vi) break;
      }
  }

 
  // Update spin configurations according to the above resolution
  for (int j = 0; j < sv.Nm; ++j) {
    if (last[j] != -1) {
      if (X[last[j]] == -1){
          lattice[j].flip();
      }
    } else {
      if (rann() < 0.5) lattice[j].flip();
    }
    last[j] = -1;
    frst[j] = -1;
  }

  delete[] lpos;
  delete[] X;
  delete[] flag;
}
