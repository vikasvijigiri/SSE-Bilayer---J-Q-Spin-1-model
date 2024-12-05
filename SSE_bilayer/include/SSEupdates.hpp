#ifndef _SSEUPDATES_HPP_DEFINED_
#define _SSEUPDATES_HPP_DEFINED_

class SSEupdates {
    	public: 
    	SSEupdates(){}

      // ********************************************************************
      // 		Update related structures 
      // ********************************************************************


      // update variables
	    int **str;
      long *frst, *last; 
     	
	    
	    const int type[3] = {1, 4, 2};
      int Nd[3], flip[5][2];
	    
     	double amax;
     	double wgt[2][2];
     	double awgt[2][2];
     	double dwgt[2][2];
	    

      // Update functions
	    void weights();


      void initialize(SSEvariables& );
	    void mcstep(SSEvariables& , SSElattice* );
	    void mcstep_measurement(SSEvariables& , SSElattice* , SSEobservables& );
	    void checkl(SSEvariables&);
	    
	    
	    void diag_update(SSEvariables& , SSElattice* );				
	    void looper(SSEvariables& , SSElattice* );
	    
      void diag_update_measurement(SSEvariables& , SSElattice* , SSEobservables& );
	    void looper_measurement(SSEvariables& , SSElattice* , SSEobservables& );
};
#endif

