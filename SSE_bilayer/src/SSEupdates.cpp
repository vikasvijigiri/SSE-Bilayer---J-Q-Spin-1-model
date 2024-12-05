
#include <iostream>

#include "../include/SSEvariables.hpp"
#include "../include/SSElattice.hpp"
#include "../include/SSEobservables.hpp"
#include "../include/SSEupdates.hpp"


static ran rann;
using namespace std;

//********************************************************************************************


void
SSEupdates::checkl(SSEvariables& sv)
{
	long Lc_new = static_cast<long>(sv.n1 + sv.n1 / 3);

	if (Lc_new > sv.Lc)
	{
		int ns = 4;
		int** tmp;
		tmp = new int *[ns];
		for (long i = 0; i < ns; ++i) tmp[i] = new int[Lc_new];
		
		// Store in a temporary array.
		for (int j = 0; j < ns; ++j)
		{
			for (long i = 0; i < Lc_new; ++i)
			{
				if (i < sv.Lc)
				{
					tmp[j][i] = str[j][i];
				}
				else
				{
					tmp[j][i] = 0;
				}
			}
		}
		// Delete old string.
		for (int i = 0; i < ns; ++i)
		{
			delete[] str[i];
		}
		delete[] str;

		
		// Allot values to newly created arrays.
		sv.Lc = Lc_new;
		str = new int *[ns];
		for (int i = 0; i < ns; ++i) str[i] = new int[sv.Lc];

		for (int j = 0; j < ns; ++j)
		{
			for (long i = 0; i < sv.Lc; ++i)
			{
				str[j][i] = tmp[j][i];
			}
		}

		for (int i = 0; i < ns; ++i)
		{
			delete[] tmp[i];
		}
		delete[] tmp;

	}
}




//********************************************************************************************
void
SSEupdates::weights()
{
	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			wgt[s1][s2] = 2.0 * 0.25 * pow(-1, 1 + s1) *pow(-1, 1 + s2);
			if (wgt[s1][s2] > amax)
				amax = wgt[s1][s2];
		}
	}

	// amax[j]+=1.0f;
	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			wgt[s1][s2] = amax - wgt[s1][s2];
		}
	}

	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			awgt[s1][s2] = wgt[s1][s2];
			if (awgt[s1][s2] > 1e-6)
			{
				dwgt[s1][s2] = 1.0 / awgt[s1][s2];
			}
			else
			{
				dwgt[s1][s2] = 1.e6;
			}
		}
	}
	

	// vertex type after flipping
	//  1 <---> 2
	flip[1][0] = 2;
	flip[1][0] = 2;	
	
	flip[2][0] = 1;	
	flip[2][0] = 1;
	
	//  1 <---> 3
	flip[1][1] = 3;	
	flip[1][1] = 3;	

	flip[3][1] = 1;	
	flip[3][1] = 1;	

	//  2 <---> 4
	flip[2][1] = 4;	
	flip[2][1] = 4;	
	
	flip[4][1] = 2;	
	flip[4][1] = 2;			

	//  3 <---> 4
	flip[4][0] = 3;	
	flip[4][0] = 3;	
	
	flip[3][0] = 4;	
	flip[3][0] = 4;		
}


void SSEupdates::initialize(SSEvariables& sv) {

    // H-bond part
    str = new int *[4]();
    for (int i = 0; i < 4; ++i) {
        str[i] = new int[sv.Lc]();
    }

    // Initialize H-bond part to zero
    for (int i = 0; i < 4; ++i) {
        std::fill(str[i], str[i] + sv.Lc, 0);
    }

    // Template S-bonds
    frst = new long[sv.Nm]();
    last = new long[sv.Nm]();

    // Initialize Template S-bonds to -1
    std::fill_n(frst, sv.Nm, -1);
    std::fill_n(last, sv.Nm, -1);

    std::fill_n(Nd, 3, 0);
}

