#include <fstream>
#include <iostream>
#include <iomanip>

#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"




// constants definition
const int H_BONDS = 4;

// Local Functions declaration
void set_H_Bonds(int** , int , int , int , int , int , int* );//, double* , double);
void set_QQ_Bonds(int** , int , int , int , int , int , int , int , int* );//, double* , double);
void set_Bi_Bonds(int** , int , int , int , int , int , int* );//, double* , double);



// Function to declare variables and allot them the values from input file.
void
SSEvariables::declare_variables()
{

  // ***************************************************************************************
  //                  Read variables from file "input_param.dat" 
  // ***************************************************************************************
	// Read input variables from file.
	std::string INPUT_FILE = "input_param.dat";
	std::ifstream vars(INPUT_FILE);
	vars >> lx;             // Size along x
	vars >> ly;             // Size along y
	vars >> JH; 			      // Heisenberg AFM interaction strength
	vars >> JQ;             // QQ AFM interaction strength
	vars >> JB; 			      // Biquad AFM interaction strength	
	vars >> Beta;           // Inverse temp
	vars >> isteps;         // Equilibration steps
	vars >> iter;           // Measurement steps
  vars >> lambda;         // Vertical coupling constant
  vars >> time_series;    // You want to print time_series data?
	
	
	Ns = lx * ly;				// No. of spins in a single layer.
	Nm = 2 * Ns;		    // No. of spins totally.

	// For PBC square lattice. For sizes < 3x3, double counting is implied.	
	NH = 8*Ns;			    // No. of Heisenberg bonds
	NQ = 8*Ns; 		      // No. of Q_Q bonds
	NB = 4*Ns;			    // No. of biquad bonds	

	Lc = 30;            // Initial SSE cutoff length
	nl = 5;          // Initial no. of loops
	n1 = 0;             // Initial no. of diagonal operators


  // ***************************************************************************************
  //                  Probabilities to insert and delete
  // ***************************************************************************************

  double coef_JH = 0.5 * JH * NH * Beta; 
  double coef_JQ = 0.5 * 0.5 * JQ * NQ * Beta;
  double coef_JB = 0.5 * JB * NB * Beta;


  

  prob_in = cum_prob[0] = cum_prob[3] = coef_JH / 4.;
  cum_prob[1] = cum_prob[2] = coef_JH * lambda / 4.;
  
  cum_prob[4] = coef_JQ /  4.;
  cum_prob[5] = cum_prob[6] =  coef_JQ * lambda / 4.;
  cum_prob[7] = coef_JQ * lambda * lambda / 4.;

  cum_prob[8] = coef_JB / 2.;
  cum_prob[9] = coef_JB * lambda / 2.;

  for (int i = 1; i < 10; ++i) {
      prob_in += cum_prob[i]; // Probability of insertion of diagonal operators.
  }

  for (int i = 0; i < 10; ++i) {
      cum_prob[i] /= prob_in;
  }
  prob_rm = 1.0 / prob_in;    // Probability of removal of diagonal operators
}




// Function to store information of stiffness, non-uniformity, bonds information into respective arrays.  
void SSEvariables::lattice_sites()
{

    // ***************************************************************************************
    //                 stiffness related Info
    // ***************************************************************************************

    JHsgnx = new int[8 * Ns];
    JHsgny = new int[8 * Ns];

    std::fill(JHsgnx, JHsgnx + 4*Ns, 1);
    std::fill(JHsgny, JHsgny + 4*Ns, 0);
    std::fill(JHsgnx + 4*Ns, JHsgnx + 8 * Ns, 0);
    std::fill(JHsgny + 4*Ns, JHsgny + 8 * Ns, 1);

    JQsgnx = new int[8 * Ns];
    JQsgny = new int[8 * Ns];

    std::fill(JQsgnx, JQsgnx + 4*Ns, 1);
    std::fill(JQsgny, JQsgny + 4*Ns, 0);
    std::fill(JQsgnx + 4*Ns, JQsgnx + 8 * Ns, 0);
    std::fill(JQsgny + 4*Ns, JQsgny + 8 * Ns, 1);

    JBsgnx = new int[4 * Ns];
    JBsgny = new int[4 * Ns];

    std::fill(JBsgnx, JBsgnx + 2*Ns, 1);
    std::fill(JBsgny, JBsgny + 2*Ns, 0);
    std::fill(JBsgnx + 2*Ns, JBsgnx + 4 * Ns, 0);
    std::fill(JBsgny + 2*Ns, JBsgny + 4 * Ns, 1);



    // ***************************************************************************************
    //                lambda_\alpha,\beta: Non uniformity of the bond coupling
    // ***************************************************************************************


//    Hsgn = new double[8 * Ns];
//    QQsgn = new double[8 * Ns];
//    Bsgn = new double[4 * Ns];


    // ***************************************************************************************
    //                Bond Info
    // ***************************************************************************************

    JHsites = new int*[8 * Ns];
    JQsites = new int*[8 * Ns];
    JBsites = new int*[4 * Ns];

    for (int i = 0; i < 8 * Ns; ++i) {
        JHsites[i] = new int[2];
        JQsites[i] = new int[8];
    }

    for (int i = 0; i < 4 * Ns; ++i) {
        JBsites[i] = new int[4];
    }

    for (int y1 = 0; y1 < ly; ++y1) {
        for (int x1 = 0; x1 < lx; ++x1) {
            // ***************************************************************************************
            //       Heisenberg interaction bonds
            // ***************************************************************************************

            // Horizontal bonds
            int s1 = x1 + y1 * lx;
            int x2 = (x1 + 1) % lx;
            int s2 = x2 + y1 * lx;

            int baseIndexHorizontal = 4 * s1;

            int shiftValues13[2] = {0, 0}; // lambda = 1 
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues13);//, Hsgn, 1.);            
            
            int shiftValues14[2] = {0, 1}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 1, Ns, shiftValues14);//, Hsgn, lambda);           

            int shiftValues15[2] = {1, 0}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 2, Ns, shiftValues15);//, Hsgn, lambda);           

            int shiftValues16[2] = {1, 1}; // lambda = 1
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 3, Ns, shiftValues16);//, Hsgn, 1.);          


            // vertical bonds
            x2 = x1;
            int y2 = (y1 + 1) % ly;
            s2 = x2 + y2 * lx;

            int baseIndexVertical = 4 * Ns + 4 * s1;

            int shiftValues17[2] = {0, 0}; // lambda = 1
            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 0, Ns, shiftValues17);//, Hsgn, 1.);              
            
            int shiftValues18[2] = {0, 1}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 1, Ns, shiftValues18);//, Hsgn, lambda);              

            int shiftValues19[2] = {1, 0}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 2, Ns, shiftValues19);//, Hsgn, lambda);            

            int shiftValues20[2] = {1, 1}; // lambda = 1
            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 3, Ns, shiftValues20);//, Hsgn, 1.);             





            // ***************************************************************************************
            //       QoQ interaction bonds
            // ***************************************************************************************

            // Horizontal bonds
            s1 = x1 + y1 * lx;
            x2 = (x1 + 1) % lx;
            s2 = x2 + y1 * lx;

            int s3 = (s1 + lx) % Ns;
            int s4 = (s2 + lx) % Ns;

            baseIndexHorizontal = 4 * s1;

            int shiftValues[8] = {0, 0, 1, 1, 0, 0, 1, 1}; // lambda = 1
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 0, Ns, shiftValues);//, QQsgn, 1.);           

            int shiftValues2[8] = {0, 0, 1, 1, 0, 1, 1, 0}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 1, Ns, shiftValues2);//, QQsgn, lambda);          

            int shiftValues3[8] = {0, 1, 1, 0, 0, 0, 1, 1}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 2, Ns, shiftValues3);//, QQsgn, lambda);          

            int shiftValues4[8] = {0, 1, 1, 0, 0, 1, 1, 0}; // lambda = 0.25
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 3, Ns, shiftValues4);//, QQsgn, lambda*lambda);          

            // Vertical bonds
            s1 = x1 + y1 * lx;
            s2 = (s1 + lx) % Ns;

            s3 = x2 + y1 * lx;
            s4 = (s3 + lx) % Ns;

            baseIndexVertical = 4 * Ns + 4 * s1;

            int shiftValues5[8] = {0, 0, 1, 1, 0, 0, 1, 1}; // lambda = 1
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 0, Ns, shiftValues5);//, QQsgn, 1.);            

            int shiftValues6[8] = {0, 0, 1, 1, 0, 1, 1, 0}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 1, Ns, shiftValues6);//, QQsgn, lambda);           

            int shiftValues7[8] = {0, 1, 1, 0, 0, 0, 1, 1}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 2, Ns, shiftValues7);//, QQsgn, lambda);           

            int shiftValues8[8] = {0, 1, 1, 0, 0, 1, 1, 0}; // lambda = 0.25  
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 3, Ns, shiftValues8);//, QQsgn, lambda*lambda);          

            // ***************************************************************************************
            //       Biquad interaction bonds
            // ***************************************************************************************

            // Horizontal bonds
            s1 = x1 + y1 * lx;
            x2 = (x1 + 1) % lx;
            s2 = x2 + y1 * lx;

            baseIndexHorizontal = 2 * s1;

            int shiftValues9[4] = {0, 0, 1, 1};     // lambda = 1
            set_Bi_Bonds(JBsites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues9);//, Bsgn, 1.);               

            int shiftValues10[4] = {0, 1, 1, 0};    // lambda = 0.5
            set_Bi_Bonds(JBsites, baseIndexHorizontal, s1, s2, 1, Ns, shiftValues10);//, Bsgn, lambda);         

            // Vertical bonds
            x2 = x1;
            y2 = (y1 + 1) % ly;
            s2 = x2 + y2 * lx;

            baseIndexVertical = 2 * Ns + 2 * s1;

            int shiftValues11[4] = {0, 0, 1, 1};    // lambda = 1
            set_Bi_Bonds(JBsites, baseIndexVertical, s1, s2, 0, Ns, shiftValues11);//, Bsgn, 1.);                 
            int shiftValues12[4] = {0, 1, 1, 0};    // lambda = 0.5
            set_Bi_Bonds(JBsites, baseIndexVertical, s1, s2, 1, Ns, shiftValues12);//, Bsgn, lambda);                 
        }
    }
}

void set_H_Bonds(int** bonds, int baseIndex, int s1, int s2, int offset, int N, int* shift) {//, double *sgn, double lmbda)  {
    bonds[baseIndex + offset][0] = s1 + N * shift[0];
    bonds[baseIndex + offset][1] = s2 + N * shift[1];
    //sgn[baseIndex + offset] = lmbda;

}

void set_QQ_Bonds(int** bonds, int baseIndex, int s1, int s2, int s3, int s4, int offset, int N, int* shift) {//, double *sgn, double lmbda) {
    bonds[baseIndex + offset][0] = s1 + N * shift[0];
    bonds[baseIndex + offset][1] = s2 + N * shift[1];
    bonds[baseIndex + offset][2] = s1 + N * shift[2];
    bonds[baseIndex + offset][3] = s2 + N * shift[3];
    bonds[baseIndex + offset][4] = s3 + N * shift[4];
    bonds[baseIndex + offset][5] = s4 + N * shift[5];
    bonds[baseIndex + offset][6] = s3 + N * shift[6];
    bonds[baseIndex + offset][7] = s4 + N * shift[7];
    //sgn[baseIndex + offset] = lmbda;
}

void set_Bi_Bonds(int** bonds, int baseIndex, int s1, int s2, int offset, int N, int* shift) {//, double *sgn, double lmbda) {
    bonds[baseIndex + offset][0] = s1 + N * shift[0];
    bonds[baseIndex + offset][1] = s2 + N * shift[1];
    bonds[baseIndex + offset][2] = s1 + N * shift[2];
    bonds[baseIndex + offset][3] = s2 + N * shift[3];
    //sgn[baseIndex + offset] = lmbda;
}


