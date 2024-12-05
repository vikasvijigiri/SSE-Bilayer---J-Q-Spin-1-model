#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>



#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"
#include "../include/SSEobservables.hpp"



using namespace std;

// Local abstract structure of the lattice sites. 
struct Site {
    int x, y;
};



// Local helper functions that helps to create position vectors.
Site* calculatePositionVectors(int , int );

// Local helper functions that helps to create dynamic arrays.
double**** create4DArray(int , int , int , int ); 
double** create2DArray(int , int );




// ****************************************************************************
//                                  Main
// ****************************************************************************



// Function that processes all the observables.
void SSEobservables::observables_processing(SSEvariables& sv) {


    Site* psite = calculatePositionVectors(sv.lx, sv.ly);
    
    // ************************** Initializer **********************************

//    // Neel order params;
//    for (auto& subarray1D : N_real) {
//        std::fill(std::begin(subarray1D), std::end(subarray1D), 0.);
//    }

//    // VBS order params; No. 1
//    for (auto& subarray2D : V_real) {
//        for (auto& subarray1D : subarray2D) {
//            std::fill(std::begin(subarray1D), std::end(subarray1D), 0.);
//        }
//    }

//    // VBS order params; No. 2
//    for (auto& subarray2D : B_real) {
//        for (auto& subarray1D : subarray2D) {
//            std::fill(std::begin(subarray1D), std::end(subarray1D), 0.);
//        }
//    }

    // VBS No.2 order params;

    // *************************** Processing ***********************************
    // Outer loop over all spins/sites
//    for (int j = 0; j < sv.Ns; ++j) {
//        double x = psite[j].x;
//        double y = psite[j].y;

//        for (int k = 0; k < 4; ++k) {
//            double V = B[k][j] / (1. * sv.iter * sv.Ns * sv.Ns); 
//            N_real[0][k] += cos(x * pi + y * (pi + 2. * pi / sv.ly)) * V;         // Along Y: x,y+dy
//            N_real[1][k] += cos(x * (pi + 2. * pi / sv.lx) + y * pi) * V;         // Along X: x+dx,y
//            N_real[2][k] += cos(x * pi + y * pi) * V;                             // Along XY: x,y
//            N_real[3][k] += cos(x * pi + y * pi) * V;                             // Along XY: x,y
//        }



//        // Inner loop over spins/sites
//        for (int i = j; i < sv.Ns; ++i) {
//            double x = psite[j].x - psite[i].x;
//            double y = psite[j].y - psite[i].y;

//            auto updateV = [&](double& O1, double Cij, double a, double b, double J) {
//                O1 += cos(x * a + y * b) * Cij / (sv.Beta * sv.Beta * J) / (1. * sv.iter * sv.Ns * sv.Ns);
//            };

//            // V_VBS;  VBS order No.1 part
//            for (int k = 0; k < 4; ++k) {      
//                for (int l = 0; l < 4; ++l) {
//                  double H_lmbda1 = (((4*i + k) % 4 == 1 || (4*i + k) % 4 == 2) && sv.lambda > 1e-6) ? sv.JH * sv.lambda : sv.JH;
//                  double H_lmbda2 = (((4*j + l) % 4 == 1 || (4*j + l) % 4 == 2) && sv.lambda > 1e-6) ? sv.JH * sv.lambda : sv.JH;
//                  double lmbda = H_lmbda1*H_lmbda2;
//                  updateV(V_real[0][k][l], Cx[4*i+k][4*j+l],  pi, 2. * pi/sv.lx, lmbda); // x1
//                  updateV(V_real[1][k][l], Cx[4*i+k][4*j+l],  pi, 0.,  lmbda);            // x2
//                  updateV(V_real[2][k][l], Cy[4*i+k][4*j+l],  2. * pi/sv.ly, pi, lmbda); // y1
//                  updateV(V_real[3][k][l], Cy[4*i+k][4*j+l],  0., pi, lmbda);            // y2
////            if ((k == 0 && l == 0) ){
////                std::cout << j << "   "<< i << " | " << Cx[4*i+k][4*j+l] / (sv.Beta * sv.Beta * lmbda) / 
////        static_cast<double>(sv.iter * sv.Ns * sv.Ns) << "  ::  "<<  V_real[0][k][l] << "  "<<  V_real[1][k][l] << "  " 
////                                    << V_real[2][k][l]  << "  " << V_real[3][k][l] << "   "<< cos(x * pi + y * 2. * pi/sv.lx) << "   " 
////                                    << cos(x * pi) << std::endl;}
//                }
//            }
//            // B_VBS;  VBS order No.2 part
//            for (int k = 0; k < 2; ++k) {      
//                for (int l = 0; l < 2; ++l) {
//                  double B_lmbda1 = (((2*i + k) % 2 == 1) && sv.lambda > 1e-6) ? sv.JB * sv.lambda : sv.JB;
//                  double B_lmbda2 = (((2*j + l) % 2 == 1) && sv.lambda > 1e-6) ? sv.JB * sv.lambda : sv.JB;
//                  double lmbda = B_lmbda1*B_lmbda2;
//                  updateV(B_real[0][k][l], Dx[2*i+k][2*j+l],  pi, 2. * pi/sv.lx, lmbda);  // x1
//                  updateV(B_real[1][k][l], Dx[2*i+k][2*j+l],  pi, 0., lmbda);             // x2
//                  updateV(B_real[2][k][l], Dy[2*i+k][2*j+l],  2. * pi/sv.ly, pi, lmbda);  // y1
//                  updateV(B_real[3][k][l], Dy[2*i+k][2*j+l],  0., pi, lmbda);             // y2
//                }
//            }

////            // B_VBS;  Symmetric VBS No. 2 order parameter. very similar to spin-1
////            updateV(B_real_pp[0], Ex[i][j],  pi, 2. * pi/sv.lx, (sv.JB*sv.JB)*(1+sv.lambda)*(1+sv.lambda));  // x1
////            updateV(B_real_pp[1], Ex[i][j],  pi, 0., sv.JB*sv.JB);             // x2
////            updateV(B_real_pp[2], Ey[i][j],  2. * pi/sv.ly, pi, sv.JB*sv.JB);  // y1
////            updateV(B_real_pp[3], Ey[i][j],  0., pi, sv.JB*sv.JB);

//        }
//    }


//    for (int i = 0; i < 4*sv.Ns; ++i) {
//        for (int j = 0; j < 4*sv.Ns; ++j) {
//                  double H_lmbda1 = ((i % 4 == 1 || i % 4 == 2) && sv.lambda > 1e-6) ? sv.JH * sv.lambda : sv.JH;
//                  double H_lmbda2 = ((j % 4 == 1 || j % 4 == 2) && sv.lambda > 1e-6) ? sv.JH * sv.lambda : sv.JH;
//                  double lmbda = H_lmbda1*H_lmbda2;

//                  std::cout << i << ", "<< j << " xx  " << "[" << sv.JHsites[i][0] << " " << sv.JHsites[i][1] << "]" << "  "
//                            << "[" << sv.JHsites[j][0] << " "  << sv.JHsites[j][1] << "]" << "   ->   "
//                            << Cx[i][j] / (sv.Beta * sv.Beta * lmbda) / static_cast<double>(sv.iter * sv.Ns * sv.Ns) << "\n";
//        }
//    }

//    for (int i = 0; i < 2*sv.Ns; ++i) {
//        for (int j = 0; j < 2*sv.Ns; ++j) {
//                  double B_lmbda1 = (i % 2 == 1) ? sv.JB * sv.lambda : sv.JB;
//                  double B_lmbda2 = (j % 2 == 1) ? sv.JB * sv.lambda : sv.JB;

//                  std::cout << i << ", "<< j << " xx  " << "[" << sv.JBsites[i][0] << " " << sv.JBsites[i][1] << "]" << "  "
//                            << "[" << sv.JBsites[i][2] << " " << sv.JBsites[i][3] << "]" << "  "
//                            << "[" << sv.JBsites[j][0] << " " << sv.JBsites[j][1] << "]" << "  "
//                            << "[" << sv.JBsites[j][2] << " " << sv.JBsites[j][3] << "]" << "   ->   "
//                            << Dx[i][j] / (sv.Beta * sv.Beta * B_lmbda1*B_lmbda2) / static_cast<double>(sv.iter * sv.Ns * sv.Ns) << "\n";
//        }
//    }

//    int bnd1 = 0;
//    int bnd2 = 0;
//    std::cout << sv.JBsites[bnd1][0] << " " << sv.JBsites[bnd1][1] << " " << sv.JBsites[bnd1][2] << " " 
//              << sv.JBsites[bnd1][3] << " " << "\n";
//    std::cout << sv.JBsites[bnd2][0] << " " << sv.JBsites[bnd2][1] << " " << sv.JBsites[bnd2][2] << " " 
//              << sv.JBsites[bnd2][3] << " " << "\n" << "\n";

//    int k, l;

//    k = 0; l = 1;
//    double d1 = (Dx[2*bnd1+k][2*bnd2+l]) / 
//                (sv.Beta * sv.Beta * sv.JB * sv.JB) / static_cast<double>(sv.iter * sv.Ns * sv.Ns);

//    k = 0; l = 1;
//    double d2 = (Dx[2*bnd1+k][2*bnd2+l]) / 
//                (sv.Beta * sv.Beta * sv.JB * sv.JB) / static_cast<double>(sv.iter * sv.Ns * sv.Ns);

//    k = 1; l = 0;
//    double d3 = (Dx[2*bnd1+k][2*bnd2+l]) / 
//                (sv.Beta * sv.Beta * sv.JB * sv.JB) / static_cast<double>(sv.iter * sv.Ns * sv.Ns);

//    k = 1; l = 1;
//    double d4 = (Dx[2*bnd1+k][2*bnd2+l]) / 
//                (sv.Beta * sv.Beta * sv.JB * sv.JB) / static_cast<double>(sv.iter * sv.Ns * sv.Ns);

//    std::cout << d1 << std::endl;

//    std::cout << d2 << std::endl;

//    std::cout << d3 << std::endl;

//    std::cout << d4 << std::endl;

//    std::cout << d1+d2+d3+d4 << std::endl;

    int klN[2] = {0, 1};
    double factorN =  (sv.iter * sv.Ns);
    for (int t = 0; t < 2; ++t) {
        int k = klN[t];

        for (int i = 0; i < sv.Ns; ++i) {
            double x = psite[i].x;
            double y = psite[i].y;

            
            double c = B[k][i] / factorN; //+ (CN[2*i+1][2*j][0] + CN[2*i][2*j+1][0] + CN[2*i+1][2*j+1][0])

            N_mx1[t] += c * cos(x * M_PI + y * (M_PI + 2.0 * M_PI / sv.lx));
            N_mx2[t] += c * cos(x * M_PI + y * M_PI);  
            N_my1[t] += c * cos(x * (M_PI + 2.0 * M_PI / sv.ly) + y * M_PI);
            N_my2[t] += c * cos(x * M_PI + y * M_PI); 
        }

        //O_N[k] = (mx2 + my2) / 2.0;
        //R_N[k] = ((fabs(mx2) > 1e-6 && fabs(my2) > 1e-6) ? (1.0 - (mx1 / mx2 + my1 / my2) / 2.0)  : 0.0);
    }


    int klV[2] = {0, 3};
    double factorV =  (sv.Beta * sv.Beta * sv.iter * sv.Ns * sv.Ns);

    for (int t = 0; t < 2; ++t) {
        int k = klV[t];
        for (int j = 0; j < sv.Ns; ++j) {
            for (int i = 0; i < sv.Ns; ++i) {
                double x = (psite[i].x - psite[j].x);
                double y = (psite[i].y - psite[j].y);

                
                double lmbda = factorV * sv.JH * sv.JH;


                double cx = Cx[4 * i + k][4 * j + k] / lmbda; //+ (CN[2*i+1][2*j][0] + CN[2*i][2*j+1][0] + CN[2*i+1][2*j+1][0])
                double cy = Cy[4 * i + k][4 * j + k] / lmbda;


                V_mx1[t] += cx * cos(x * M_PI + y * (2.0 * M_PI / sv.lx));
                V_mx2[t] += cx * cos(x * M_PI);
                V_my1[t] += cy * cos(x * (2.0 * M_PI / sv.ly) + y * M_PI);
                V_my2[t] += cy * cos(y * M_PI);
            }
//                std::cout  << k << "  " << 4*i + k << "  |  "<< cx << "  |  " << mx1 << " " << mx2 << " "  << std::endl;
        }
        //std::cout << V_mx1[t] / sv.Ns << "  " << V_mx2[t] / sv.Ns << "  " << V_my1[t] / sv.Ns << "  " << V_my2[t] / sv.Ns << std::endl;
        //O_V[k][k] = (mx2 + my2) / 2.0;
        //R_V[k][k] = ((fabs(mx2) > 1e-6 && fabs(my2) > 1e-6) ? (1.0 - (mx1 / mx2 + my1 / my2) / 2.0)  : 0.0);
    }

    int klB[2] = {0, 1};
    double JB_factor[2] = {sv.JB, sv.JB*sv.lambda};
    double factorB =  (sv.Beta * sv.Beta * sv.iter * sv.Ns * sv.Ns);
    
    for (int t = 0; t < 2; ++t) {
        int k = klB[t];
        for (int j = 0; j < sv.Ns; ++j) {
            for (int i = 0; i < sv.Ns; ++i) {
                double x = (psite[i].x - psite[j].x);
                double y = (psite[i].y - psite[j].y);


                double lmbda = factorB * JB_factor[k] * JB_factor[k];

                double dx = Dx[2 * i + k][2 * j + k] / lmbda; //+ (CN[2*i+1][2*j][0] + CN[2*i][2*j+1][0] + CN[2*i+1][2*j+1][0])
                double dy = Dy[2 * i + k][2 * j + k] / lmbda;

                B_mx1[t] += dx * cos(x * M_PI + y * (2.0 * M_PI / sv.lx));
                B_mx2[t] += dx * cos(x * M_PI);
                B_my1[t] += dy * cos(x * (2.0 * M_PI / sv.ly) + y * M_PI);
                B_my2[t] += dy * cos(y * M_PI);
            }
        }
        //O_B[k][k] = (mx2 + my2) / 2.0;
        //R_B[k][k] = ((fabs(mx2) > 1e-6 && fabs(my2) > 1e-6) ? (1.0 - (mx1 / mx2 + my1 / my2) / 2.0)  : 0.0);
    }


//    double B_realx1_sym = B_real[0][0][0] + (B_real[0][0][1] + B_real[0][1][0] + B_real[0][1][1]);
//    double B_realx2_sym = B_real[1][0][0] + (B_real[1][0][1] + B_real[1][1][0] + B_real[1][1][1]);
//    double B_realy1_sym = B_real[2][0][0] + (B_real[2][0][1] + B_real[2][1][0] + B_real[2][1][1]);
//    double B_realy2_sym = B_real[3][0][0] + (B_real[3][0][1] + B_real[3][1][0] + B_real[3][1][1]);
    // Post processing part for calculating order parameter O_N and R_N. For this case, it is done later in python.
//    auto processOrderParams = [&](double& O,
//                                  double& R, double& Real1, double& Real2, double& Real3, double& Real4) {
//        double O1 = Real1;
//        double O2 = Real2;
//        double O3 = Real3; 
//        double O4 = Real4;
//        O = (O3 + O4) / 2.;
//        R = 1. - (O1 / O3 + O2 / O4) / 2.;
//    };

//    // V_VBS;  VBS order No.1 part
//    for (int k = 0; k < 4; ++k) {      
//        double realx1 = N_real[0][k];
//        double realy1 = N_real[1][k];
//        double realx2 = N_real[2][k];
//        double realy2 = N_real[3][k];
//        processOrderParams(O_N[k], R_N[k], realx1, realy1, realx2, realy2);
//    }

//    // V_VBS;  VBS order No.1 part
//    for (int k = 0; k < 4; ++k) {      
//        for (int l = 0; l < 4; ++l) {
//            double realx1 = V_real[0][k][l];
//            double realx2 = V_real[1][k][l];
//            double realy1 = V_real[2][k][l];
//            double realy2 = V_real[3][k][l];
//            processOrderParams(O_V[k][l], R_V[k][l], realx1, realy1, realx2, realy2);
//            if ((k == 0 && l == 0) || (k == 3 && l == 3))
//                std::cout << k << "   "<< l << " | "<<  realx1 << "  "<<  realy1 << "  " << realx2  << "  " << realy2 << std::endl;
//        }
//    }

//    // V_VBS;  VBS order No.1 part
//    for (int k = 0; k < 2; ++k) {      
//        for (int l = 0; l < 2; ++l) {
//            double realx1 = B_real[0][k][l];
//            double realx2 = B_real[1][k][l];
//            double realy1 = B_real[2][k][l];
//            double realy2 = B_real[3][k][l];
//            processOrderParams(O_B[k][l], R_B[k][l], realx1, realy1, realx2, realy2);
//        }
//    }


//    // form :                   (= =)      +        (= x)         +     (x =)       +    (x x)
//    double B_realx1_asym = B_real[0][0][0] + -1.*(B_real[0][0][1] + B_real[0][1][0] + B_real[0][1][1]);
//    double B_realx2_asym = B_real[1][0][0] + -1.*(B_real[1][0][1] + B_real[1][1][0] + B_real[1][1][1]);
//    double B_realy1_asym = B_real[2][0][0] + -1.*(B_real[2][0][1] + B_real[2][1][0] + B_real[2][1][1]);
//    double B_realy2_asym = B_real[3][0][0] + -1.*(B_real[3][0][1] + B_real[3][1][0] + B_real[3][1][1]);
//    //double O_B, R_B;
//    //processOrderParams(O_N1, O_N2, O_N3, O_N3, O_N, R_N, N_real1, N_real2, N_real3, N_real3);
//    //processOrderParams(O_V1, O_V2, O_V3, O_V4, O_V, R_V, V_realy1, V_realx1, V_realy2, V_realx2);

//    processOrderParams(O_B1, R_B1, B_real[2][0][0], B_real[0][0][0], B_real[3][0][0], B_real[1][0][0]);
//    processOrderParams(O_B2, R_B2, B_real[2][0][1], B_real[0][0][1], B_real[3][0][1], B_real[1][0][1]);
//    processOrderParams(O_B3, R_B3, B_real[2][1][0], B_real[0][1][0], B_real[3][1][0], B_real[1][1][0]);
//    processOrderParams(O_B4, R_B4, B_real[2][1][1], B_real[0][1][1], B_real[3][1][1], B_real[1][1][1]);    

//    processOrderParams(O_B_sym, R_B_sym, B_realy1_sym, B_realx1_sym, B_realy2_sym, B_realx2_sym);
//    processOrderParams(O_B_asym, R_B_asym, B_realy1_asym, B_realx1_asym, B_realy2_asym, B_realx2_asym);


    //std::cout << O_B_sym << "  split   " << R_B_sym << std::endl;

//    // Pre-processed as before
//    processOrderParams(O_B_sym, R_B_sym, B_real_pp[2], B_real_pp[0], B_real_pp[3], B_real_pp[1]);


//    std::cout << O_B_sym << "  No split   " << R_B_sym << std::endl;
    //std::cout << B_realx1 << " " << B_realx2 << " " << B_realy1 << " " << B_realy2 << " | " << O_B << "   " << R_B << "   " << std::endl;

    // <Energy>, <Energy^2>, <Energy^4> 
    enrg /= sv.iter;
    enrg2 /= sv.iter;
    enrg4 /= sv.iter;    

    // <Mag> 
    SMagL1 /= sv.iter; // Layer 1
    SMagL2 /= sv.iter; // Layer 2
    
    // <Mag^2>
    SMag_squareL1 /= sv.iter; // Layer 1
    SMag_squareL2 /= sv.iter; // Layer 2
    SMag_squareI1 /= sv.iter; // Inter layer
  
    // <Mag^4>
    SMag_fourL1 /= sv.iter; // Layer 1
    SMag_fourL2 /= sv.iter; // Layer 2
    SMag_fourI1 /= sv.iter; // Inter layer
    SMag_fourI2 /= sv.iter; // Inter layer
    SMag_fourI3 /= sv.iter; // Inter layer
    
    // <rho_x>, <rho_y>
    stiffx /= sv.iter;  // stiffness along x-direction
    stiffy /= sv.iter;  // stiffness along x-direction

}


// Function to create a dynamic 4D array and initialize to zero
double**** create4DArray(int size1, int size2, int size3, int size4) {
    double**** array = new double***[size1];
    for (int i = 0; i < size1; ++i) {
        array[i] = new double**[size2];
        for (int j = 0; j < size2; ++j) {
            array[i][j] = new double*[size3];
            for (int k = 0; k < size3; ++k) {
                array[i][j][k] = new double[size4]();
            }
        }
    }
    return array;
}

// Function to create a dynamic 2D array and initialize to zero
double** create2DArray(int size1, int size2) {
    double** array = new double*[size1];
    for (int i = 0; i < size1; ++i) {
        array[i] = new double[size2]();
    }
    return array;
}

// Function that initializes all the observables to zero (for measurement). 
void SSEobservables::Initiate_observables(SSEvariables& sv) {
    enrg = 0.;
    enrg2 = 0.;
    enrg4 = 0.;

    // Magnetic fields
    SMagL1 = 0.;
    SMagL2 = 0.;
    SMag_squareL1 = 0.;
    SMag_squareL2 = 0.;
    SMag_squareI1 = 0.;

    // More magnetic fields
    SMag_fourL1 = 0.;
    SMag_fourL2 = 0.;
    SMag_fourI1 = 0.;
    SMag_fourI2 = 0.;
    SMag_fourI3 = 0.;

    // Stiffness
    stiffx = 0.;
    stiffy = 0.;

    
    for (int t=0; t<2; ++t){
        // Neel params
        N_mx1[t] = 0.;
        N_mx2[t] = 0.;
        N_my1[t] = 0.;
        N_my2[t] = 0.;


        // VBS  No. 1 params
        V_mx1[t] = 0.;
        V_mx2[t] = 0.;
        V_my1[t] = 0.;
        V_my2[t] = 0.;


        // VBS  No. 2 params
        B_mx1[t] = 0.;
        B_mx2[t] = 0.;
        B_my1[t] = 0.;
        B_my2[t] = 0.;
    }

//    // Define 4D vectors (WARNING: VECTORS ARE NOT WORKING IN CHANDRA CLUSTER (IITB))
//    Cx.resize(4, std::vector<std::vector<std::vector<double>>>(4, std::vector<std::vector<double>>(sv.Ns, std::vector<double>(sv.Ns, 0.0))));
//    Cy.resize(4, std::vector<std::vector<std::vector<double>>>(4, std::vector<std::vector<double>>(sv.Ns, std::vector<double>(sv.Ns, 0.0))));
//    Dx.resize(2, std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(sv.Ns, std::vector<double>(sv.Ns, 0.0))));
//    Dy.resize(2, std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(sv.Ns, std::vector<double>(sv.Ns, 0.0))));


//    B.resize(4, std::vector<double>(sv.Ns, 0.0));

    // Convert to dynamic arrays
    Cx = create2DArray(4*sv.Ns, 4*sv.Ns);
    Cy = create2DArray(4*sv.Ns, 4*sv.Ns);
    Dx = create2DArray(2*sv.Ns, 2*sv.Ns);
    Dy = create2DArray(2*sv.Ns, 2*sv.Ns);
//    Ex = create2DArray(2*sv.Ns, 2*sv.Ns);
//    Ey = create2DArray(2*sv.Ns, 2*sv.Ns);
    B = create2DArray(4, sv.Ns);

}


// Function to calculate position vectors with periodic boundary conditions
Site* calculatePositionVectors(int Lx, int Ly) {
    Site* positions = new Site[Lx * Ly];

//    // Reference position for site 0
//    int refX = 0, refY = 0;

//    // Calculate position vectors
//    for (int i = 0; i < Lx * Ly; ++i) {
//        int relX = i % Lx;
//        int relY = i / Lx;
//        positions[i].x = relX - refX;
//        positions[i].y = relY - refY;

//        // Apply periodic boundary conditions along x
//        if (positions[i].x > Lx / 2) {
//            positions[i].x -= Lx;
//        }
//        if (positions[i].x < -Lx / 2) {
//            positions[i].x += Lx;
//        }

//        // Apply periodic boundary conditions along y
//        if (positions[i].y > Ly / 2) {
//            positions[i].y -= Ly;
//        }
//        if (positions[i].y < -Ly / 2) {
//            positions[i].y += Ly;
//        }
//    }

    int z1=0, y1=0;
	  for (int k = 1; k <= Lx*Ly; k++)
	  {
		  //Block 1(Bottom left)
		  if (k % Lx <= int((Lx) / 2) && k <= int((Lx*Ly / 2) + Lx))
		  {
			  if (k % Lx != 0)
			  {
				  z1 = k % Lx - 1;
				  y1 = int((k - 1) / Lx);
				  //g1 = z1 + y1;
				  //g2 = y1 - z1;
			  }
			  else
			  {
				  z1 = -k % Lx - 1;
				  y1 = int((k - 1) / Lx);
				  //g1 = z1 + y1;
				  //g2 = y1 - z1;
			  }
		  }
		  //Block 2(Bottom right)
		  else if (k % Lx > int((Lx) / 2) && k <= int((Lx*Ly / 2) + Lx))
		  {
			  z1 = k % Lx - 1 - Lx;
			  y1 = int((k - 1) / Lx);
			  //g1 = z1 + y1;
			  //g2 = y1 - z1;
		  }
		  //Block 3(top left)
		  else if (k % Lx <= int((Lx) / 2) && k >= int((Lx*Ly / 2) + Lx))
		  {
			  z1 = k % Lx - 1;
			  y1 = int((k - 1) / Lx) - Lx;
			  //g1 = z1 + y1;
			  //g2 = y1 - z1;
		  }
		  //Block 4(top right)
		  else if (k % Lx > int((Lx) / 2) && k >= int((Lx*Ly / 2) + Lx))
		  {
			  z1 = k%Lx - 1 - Lx;
			  y1 = int((k - 1) / Lx) - Lx;
			  //g1 = z1 + y1;
			  //g2 = y1 - z1;
		  }
		  positions[k - 1].x = z1;//g1 - 1;
		  positions[k - 1].y = y1;//g2 - 1;
	  }



    return positions;
}



const char** SSEobservables::createHeadersArray() {
    static const char* headers[33];

    // energy headers
    headers[0] = "enrg";
    headers[1] = "enrg2";
    headers[2] = "enrg4";

//    // Staggered Magnetization headers
//    headers[3] = "SMagL1";
//    headers[4] = "SMagL2";

    // Staggered Magnetization^2 headers
    headers[3] = "SMag_square_L1";
    headers[4] = "SMag_square_L2";
//    headers[5] = "SMag_squareI1";

    // Staggered Magnetization^4 headers
    headers[5] = "SMag_four_L1";
    headers[6] = "SMag_four_L2";
//    headers[8] = "SMag_fourI1";
//    headers[9] = "SMag_fourI2";
//    headers[10] = "SMag_fourI3";


    // Stiffness headers
    headers[7] = "stiffx";
    headers[8] = "stiffy";



    // First layer
    // Neel 
    headers[9]  = "N_mx1_0";
    headers[10] = "N_mx2_0";
    headers[11] = "N_my1_0";
    headers[12] = "N_my2_0";

    // VBS No. 1 
    headers[13] = "V_mx1_0";
    headers[14] = "V_mx2_0";
    headers[15] = "V_my1_0";
    headers[16] = "V_my2_0";

    // VBS No. 2
    headers[17] = "B_mx1_0";
    headers[18] = "B_mx2_0";
    headers[19] = "B_my1_0";
    headers[20] = "B_my2_0";


    // Second layer 
    // Neel 
    headers[21] = "N_mx1_1";
    headers[22] = "N_mx2_1";
    headers[23] = "N_my1_1";
    headers[24] = "N_my2_1";

    // VBS No. 1
    headers[25] = "V_mx1_1";
    headers[26] = "V_mx2_1";
    headers[27] = "V_my1_1";
    headers[28] = "V_my2_1";

    // VBS No. 2
    headers[29] = "B_mx1_1";
    headers[30] = "B_mx2_1";
    headers[31] = "B_my1_1";
    headers[32] = "B_my2_1";



//    // Neel 
//    headers[15] = "O_N_0";
//    headers[16] = "O_N_1";
//    headers[17] = "O_N_2";
//    headers[18] = "O_N_3";
//    headers[19] = "R_N_0";
//    headers[20] = "R_N_1";
//    headers[21] = "R_N_2";
//    headers[22] = "R_N_3";

//    // VBS No. 1
//    headers[23] = "O_V_00";
//    headers[24] = "R_V_00";
//    headers[25] = "O_V_01";
//    headers[26] = "R_V_01";
//    headers[27] = "O_V_10";
//    headers[28] = "R_V_10";
//    headers[29] = "O_V_11";
//    headers[30] = "R_V_11";
//    headers[31] = "O_V_20";
//    headers[32] = "R_V_20";
//    headers[33] = "O_V_21";
//    headers[34] = "R_V_21";
//    headers[35] = "O_V_22";
//    headers[36] = "R_V_22";
//    headers[37] = "O_V_23";
//    headers[38] = "R_V_23";
//    headers[39] = "O_V_30";
//    headers[40] = "R_V_30";
//    headers[41] = "O_V_31";
//    headers[42] = "R_V_31";
//    headers[43] = "O_V_32";
//    headers[44] = "R_V_32";
//    headers[45] = "O_V_33";
//    headers[46] = "R_V_33";

//    // VBS No. 2  
//    headers[47] = "O_B_00";
//    headers[48] = "R_B_00";
//    headers[49] = "O_B_01";
//    headers[50] = "R_B_01";
//    headers[51] = "O_B_10";
//    headers[52] = "R_B_10";
//    headers[53] = "O_B_11";
//    headers[54] = "R_B_11";



    return headers;
}



const double* SSEobservables::createValuesArray() {
    static double* values = new double[33];

    // Energy
    values[0] = enrg;
    values[1] = enrg2;
    values[2] = enrg4;

//    // Staggered Magnetization
//    values[3] = SMagL1;
//    values[4] = SMagL2;

    // Staggered Magnetization^2
    values[3] = SMag_squareL1;
    values[4] = SMag_squareL2;
//    values[5] = SMag_squareI1;

    // Staggered Magnetization^4
    values[5] = SMag_fourL1;
    values[6] = SMag_fourL2;
//    values[8] = SMag_fourI1;
//    values[9] = SMag_fourI2;
//    values[10] = SMag_fourI3;

    // Stiffness
    values[7] = stiffx;
    values[8] = stiffy;

    // First layer
    // Neel 
    values[9]  = N_mx1[0];
    values[10] = N_mx2[0];
    values[11] = N_my1[0];
    values[12] = N_my2[0];

    // VBS No. 1
    values[13] = V_mx1[0];
    values[14] = V_mx2[0];
    values[15] = V_my1[0];
    values[16] = V_my2[0];

    // VBS No. 2
    values[17] = B_mx1[0];
    values[18] = B_mx2[0];
    values[19] = B_my1[0];
    values[20] = B_my2[0];


    // Second layer
    // Neel 
    values[21] = N_mx1[1];
    values[22] = N_mx2[1];
    values[23] = N_my1[1];
    values[24] = N_my2[1];

    // VBS No. 1
    values[25] = V_mx1[1];
    values[26] = V_mx2[1];
    values[27] = V_my1[1];
    values[28] = V_my2[1];

    // VBS No. 2
    values[29] = B_mx1[1];
    values[30] = B_mx2[1];
    values[31] = B_my1[1];
    values[32] = B_my2[1];


//    // Neel 
//    values[15] = O_N[0];
//    values[16] = O_N[1];
//    values[17] = O_N[2];
//    values[18] = O_N[3];
//    values[19] = R_N[0];
//    values[20] = R_N[1];
//    values[21] = R_N[2];
//    values[22] = R_N[3];

//    // VBS No. 1
//    values[23] = O_V[0][0];
//    values[24] = R_V[0][0];
//    values[25] = O_V[0][1];
//    values[26] = R_V[0][1];
//    values[27] = O_V[1][0];
//    values[28] = R_V[1][0];
//    values[29] = O_V[1][1];
//    values[30] = R_V[1][1];
//    values[31] = O_V[2][0];
//    values[32] = R_V[2][0];
//    values[33] = O_V[2][1];
//    values[34] = R_V[2][1];
//    values[35] = O_V[2][2];
//    values[36] = R_V[2][2];
//    values[37] = O_V[2][3];
//    values[38] = R_V[2][3];
//    values[39] = O_V[3][0];
//    values[40] = R_V[3][0];
//    values[41] = O_V[3][1];
//    values[42] = R_V[3][1];
//    values[43] = O_V[3][2];
//    values[44] = R_V[3][2];
//    values[45] = O_V[3][3];
//    values[46] = R_V[3][3];

//    // VBS No. 2
//    values[47] = O_B[0][0];
//    values[48] = R_B[0][0];
//    values[49] = O_B[0][1];
//    values[50] = R_B[0][1];
//    values[51] = O_B[1][0];
//    values[52] = R_B[1][0];
//    values[53] = O_B[1][1];
//    values[54] = R_B[1][1];


    return values;
}
