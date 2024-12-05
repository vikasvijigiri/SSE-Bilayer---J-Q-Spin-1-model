#include "SSEvariables.hpp"
#include "SSElattice.hpp"
#include "SSEupdates.hpp"
#include <iostream>

static ran rann;
using namespace std;


void
SSEupdates::diag_update(int w, SSElattice *lattice)
{
	int b, k, o, ss1, ss2;
	double p, r, cp;


	/*
	//cout << "came here " << endl;
	// H/S-bonds operators
	for (int i = 0; i < Lc; ++i)
	{
		o = str[0][i]+str[1][i]; 
		//cout << "her e1 " << endl;
		//cout << "O is  " << o << "   " << str[0][i] << "   " << str[1][i] << "   " << str[2][i]  << endl;
		if (o == 0)
		{	
			r = rann();
			cp = 0.0;
			for (k = 0; k<2; k++){
				//cp += cum_prob[k];
				//if (cp > r and cum_prob[k] < 1e-6){	// J or Q?
					b = int(rann() * Nb[k]);
					p = 1.0;				
					if (k == 0){
						ss1 = (1+lattice[JHsites[b][0]].S())/2;
						ss2 = (1+lattice[JHsites[b][1]].S())/2;
						p = p * awgt[ss1][ss2];
						p = p * Beta* Nb[0] * J[0] * 0.5  / float(Lc-n1);
						//cout << "Weight  a " << p << "  " << i << "  " << b << "  " << p* prob_in / float(Lc-n1) << endl;
					}
					if (k == 1){
						for (int j = 0; j < 4; ++j){
							b = int(rann() *Nb[k]);
							p = 1.0;						
							ss1 = (1 + lattice[JQsites[b][2*j]].S()) / 2;
							ss2 = (1 + lattice[JQsites[b][2*j+1]].S()) / 2;
							p = p * awgt[ss1][ss2];
						}
						p = p * Beta *Nb[1] * J[1] * 0.5 * 0.5 / float(Lc-n1);
						//cout << "Weight  b " << p << "  " << i << "  " << b << "  " << p* prob_in / float(Lc-n1) << endl;
					}					
					//p = p * prob_in / float(Lc-n1);
					//cout << "prob is " << p << "   " << k << endl;					
					if (rann() <= p) // Bond
					{
						str[0][i] = 1; 		// Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to QQ and Bi) 
						str[1][i] = 1;			// Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to Bi)) 
						str[2][i] = k;			// Which type of operator, I, H, QQ, Biquad
						str[3][i] = b + 1; 		// Bond-index
						Nd[k] += 1;
						n1 += 1;
						break; 
					}
				//break;					
				//}	
			}
			//cout << "her e1 " << endl;		 								
		}
		else if (o == 2)
		{
			k = str[2][i];
			if (k == 0){
				p = float(Lc - n1 + 1) / ( Beta* Nb[0] * J[0] * 0.5 );
			}
			if (k == 1){
				p = float(Lc - n1 + 1) / ( Beta *Nb[1] * J[1] * 0.5 * 0.5 );
			}
			//p = prob_rm * float(Lc - n1 + 1); 
			if (rann() <= p)
			{
				Nd[k] -= 1;
				str[0][i] = 0;
				str[1][i] = 0;
				str[2][i] = 0;
				str[3][i] = 0;
				n1 -= 1;
			}
		}
		else
		{
			b = str[3][i]-1;
			k = str[2][i];
			if (k == 0){
				lattice[JHsites[b][0]].flip();
				lattice[JHsites[b][1]].flip();			
			}
			else if (k == 1){
				for (int j = 0; j < 2; j++){
					o = str[j][i];			
					if (o == 2){
						lattice[JQsites[b][0+4*j]].flip();
						lattice[JQsites[b][1+4*j]].flip();
					}
					else if (o == 3){
						lattice[JQsites[b][2+4*j]].flip();
						lattice[JQsites[b][3+4*j]].flip();
					}
					else if (o == 4){
						lattice[JQsites[b][0+4*j]].flip();
						lattice[JQsites[b][1+4*j]].flip();
						lattice[JQsites[b][2+4*j]].flip();
						lattice[JQsites[b][3+4*j]].flip();
					}
				}			
			}			
		}
		//cout << " i is    " <<  i << "  " << str[0][i] << "  " << str[1][i] << "   type is   " << str[2][i] << "  " << int(n1 + n1 / 3) << "  " << Lc << endl;
	}
	//cout << "  " << endl;
	*/
	SSElattice * lattice_temp;
	lattice_temp = new SSElattice [Nm];
	
	for (int i = 0; i < Nm; ++i)
	{
		lattice_temp[i].set_S(lattice[i].S());
	}
	
	// Seperate loop to insert the diagonal J.
	for (int i = 0; i < Lc; ++i)
	{
		o = str[0][i];
		k = str[2][i];
		if (o == 0)
		{
			b = int(rann() *(Nb[0]));
			ss1 = (1 + lattice_temp[JHsites[b][0]].S()) / 2;
			ss2 = (1 + lattice_temp[JHsites[b][1]].S()) / 2;
			p = awgt[ss1][ss2] * (0.5 * Beta * J[0] * float(Nb[0]) / float(Lc - n1));
			//if (w==10)cout << "Inserted at i = " << i << " at bond = " << b << " with spins =  "<< lattice_temp[JHsites[b][0]].S() << ","<<lattice_temp[JHsites[b][1]].S() << endl;
			if (p > 1.0 or rann() <= p)
			{
				str[0][i] = 1;
				str[1][i] = 1;
				str[2][i] = 0;
				str[3][i] = b + 1;
				n1 += 1;
				Nd[0] += 1;
			}
		}
		else if (o == 1)
		{
			if (k == 0){		
				p =  (float(Lc - n1 + 1) / (0.5 * Beta * J[0] * float(Nb[0])));
				if (p > 1.0 or rann() <= p)
				{
					Nd[k] -= 1;
					str[0][i] = 0;
					str[1][i] = 0;
					str[2][i] = 0;
					str[3][i] = 0;
					n1 -= 1;
				}
			}
		}
		else
		{
			b = str[3][i] - 1;
			if (k == 0){
				lattice_temp[JHsites[b][0]].flip();
				lattice_temp[JHsites[b][1]].flip();
			}
			else if (k == 1){
				for (int j = 0; j < 2; j++){
					o = str[j][i];			
					if (o == 2){
						lattice_temp[JQsites[b][0+4*j]].flip();
						lattice_temp[JQsites[b][1+4*j]].flip();
					}
					else if (o == 3){
						lattice_temp[JQsites[b][2+4*j]].flip();
						lattice_temp[JQsites[b][3+4*j]].flip();
					}
					else if (o == 4){
						lattice_temp[JQsites[b][0+4*j]].flip();
						lattice_temp[JQsites[b][1+4*j]].flip();
						lattice_temp[JQsites[b][2+4*j]].flip();
						lattice_temp[JQsites[b][3+4*j]].flip();
					}
				}
			}
		}
	}
	for (int i = 0; i < Np; ++i)
	{
		o = vxoper[pstr[i]];
		if (o == 2)
		{
			//b = pstr[1][i] - 1;
			//s1 = i; //Jsites[b][0];
			//s2 = i+Np; //Jsites[b][1];
			lattice[i].flip();
			lattice[i+Np].flip();
		}
	}

	// ****************************************************************************
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//  				Second loop starts here 
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	// ****************************************************************************	
	
	for (int i = 0; i < Nm; ++i)
	{
		lattice_temp[i].set_S(lattice[i].S());
	}	
	
	// Seperate loop to insert the diagonal QQ.
	for (int i = 0; i < Lc; ++i)
	{	
		// lattice_temp = lattice
		o = str[0][i]*str[1][i];
		k = str[2][i];
		if (o == 0)
		{
				b = int(rann() *(Nb[1]));
				p = 1.0;
				for (int j = 0; j < 4; ++j){		
					ss1 = (1 + lattice_temp[JQsites[b][2*j]].S()) / 2;
					ss2 = (1 + lattice_temp[JQsites[b][2*j+1]].S()) / 2;				
					p = p*awgt[ss1][ss2];
				}				
				p = p * (Beta *Nb[1] * J[1] * 0.5 * 0.5) / float(Lc-n1);
				if (rann() <= p)
				{
					str[0][i] = 1;
					str[1][i] = 1;
					str[2][i] = 1;
					str[3][i] = b + 1;
					n1 += 1;
					Nd[1] += 1;
				}
		}
		else if (o == 1)
		{
			if (k == 1){
				p =  float(Lc - n1 + 1)/(Beta *Nb[1] * J[1] * 0.5 * 0.5);
				if (rann() <= p)
				{
					Nd[k] -= 1;
					str[0][i] = 0;
					str[1][i] = 0;
					str[2][i] = 0;
					str[3][i] = 0;				
					n1 -= 1;
				}
			}	
		}
		else
		{
			b = str[3][i] - 1;
			if (k == 0){			
				lattice_temp[JHsites[b][0]].flip();
				lattice_temp[JHsites[b][1]].flip();
			}
			if (k == 1){
				for (int j = 0; j < 2; j++){
					o = str[j][i];			
					if (o == 2){
						lattice_temp[JQsites[b][0+4*j]].flip();
						lattice_temp[JQsites[b][1+4*j]].flip();
					}
					else if (o == 3){
						lattice_temp[JQsites[b][2+4*j]].flip();
						lattice_temp[JQsites[b][3+4*j]].flip();
					}
					else if (o == 4){
						lattice_temp[JQsites[b][0+4*j]].flip();
						lattice_temp[JQsites[b][1+4*j]].flip();
						lattice_temp[JQsites[b][2+4*j]].flip();
						lattice_temp[JQsites[b][3+4*j]].flip();
					}
				}
			}	
		}
	}
  for (int i = 0; i < Np; ++i)
	{
		o = vxoper[pstr[i]];
		if (o == 2)
		{
			//b = pstr[1][i] - 1;
			//s1 = i; //Jsites[b][0];
			//s2 = i+Np; //Jsites[b][1];
			lattice[i].flip();
			lattice[i+Np].flip();
		}
	}	
}

