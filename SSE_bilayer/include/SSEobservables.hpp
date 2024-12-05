#ifndef _SSEOBSERVABLES_HPP_DEFINED_
#define _SSEOBSERVABLES_HPP_DEFINED_

#pragma once
class SSEobservables {
  	public:

			  // Observables
//        double enrg, enrg2, enrg4, mag1, mag2, mag_square1, mag_square2, mag_square3;
//        double mag_four1, mag_four2, mag_four3, mag_four4, mag_four5;
//        double stiffx, stiffy;

        double enrg;
        double enrg2;
        double enrg4;

        // Magnetic fields
        double SMagL1;
        double SMagL2;
        double SMag_squareL1;
        double SMag_squareL2;
        double SMag_squareI1;

        // More magnetic fields
        double SMag_fourL1;
        double SMag_fourL2;
        double SMag_fourI1;
        double SMag_fourI2;
        double SMag_fourI3;

        // Stiffness
        double stiffx;
        double stiffy;

        double O_N[4];
        double R_N[4];

        double O_V[4][4];
        double R_V[4][4];

        double O_B[2][2];
        double R_B[2][2];

//        double O_B_sym;
//        double R_B_sym;

//        double O_B_asym;
//        double R_B_asym;

        double N_mx1[2];
        double N_mx2[2];
        double N_my1[2];
        double N_my2[2];

        double V_mx1[2];
        double V_mx2[2];
        double V_my1[2];
        double V_my2[2];

        double B_mx1[2];
        double B_mx2[2];
        double B_my1[2];
        double B_my2[2];



//        // Neel order params;
//        double N_real[4][4];

//        // VBS order params; No. 1
//        double V_real[4][4][4];


//        // VBS order params; No. 2
//        double B_real[4][2][2];



        double **B;
        double **Cx, **Cy;
        double **Dx, **Dy;
        double **Ex, **Ey; 

        // Observable functions
        void observables_processing(SSEvariables& );
        void Initiate_observables(SSEvariables& );
        const char** createHeadersArray();
        const double* createValuesArray();


};
#endif
