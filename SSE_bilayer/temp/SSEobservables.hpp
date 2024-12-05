#ifndef _SSEOBSERVABLES_HPP_DEFINED_
#define _SSEOBSERVABLES_HPP_DEFINED_
#pragma once
class SSEobservables {
  	public:
  	
  			// Variables
      	double enrg, enrg_bin1, enrg_bin2, enrg_avg, enrg_err; //*enrg_bin, 
      	double R_N, R_N_bin1, R_N_bin2, R_N_avg, R_N_err; //*R_N_bin, 
      	double O_N, O_N_bin1, O_N_bin2, O_N_avg, O_N_err; //*O_N_bin, 
      	double R_V, R_V_bin1, R_V_bin2, R_V_avg, R_V_err; //*R_V_bin, 
      	double O_V, O_V_bin1, O_V_bin2, O_V_avg, O_V_err; //*O_V_bin,      	
      	int **psite;
      	
      	
      	
      	// Functions
      	void observables(SSElattice* , SSEupdates& );
      	void Jackniffe_data(int , int , int );
      	std::pair<double, double>  Jackniffe(double *, double , int );
      	void Initiate_observables();
};
#endif
