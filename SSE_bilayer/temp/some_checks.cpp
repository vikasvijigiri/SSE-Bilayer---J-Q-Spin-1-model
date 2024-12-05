#include "../include/SSEvariables.hpp"
#include "../include/SSElattice.hpp"
#include "../include/SSEupdates.hpp"
#include <iostream>

static ran rann;
using namespace std;

void SSEupdates::print_state_string(SSElattice *lattice)
{
	SSElattice * lattice_temp;
	lattice_temp = new SSElattice [Nm];
	
	cout<<"ms config entering     : ";
	for (int i = 0; i < Nm; ++i) {
		lattice_temp[i].set_S(lattice[i].S());
		cout<<lattice_temp[i].S()<<"  ";		
		} 
	cout<<endl;

	
	//cout<<"0: iden, 1: diag, >=2: off-diag"<<endl;
	//cout<<"0: J-type, 1: QoQ"<<endl;
	//cout<<""<< endl;
	//printing state propagration and op. str. info 
	for (int i = 0; i < Lc; ++i) {

		int b = str[3][i] - 1;
    int k = str[2][i];
    if (str[0][i]*str[1][i] >= 1){
		  if (k==0){
			cout << "i: " << i <<" | I/D/OD: "<<str[0][i]<< " | J/QoQ: "<<str[2][i] <<" | sites: ("<<JHsites[b][0]<<","<<JHsites[b][1]<<") "; 
			} 
			else if (k == 1){ 
			cout << "i: " << i <<" | I/D/OD: "<<str[0][i]<<","<<str[1][i]<< " | J/QoQ: "<<str[2][i] <<" | sites: ("<<JQsites[b][0]<<","<<JQsites[b][1]<<"), "
			<< "("<<JQsites[b][2]<<","<<JQsites[b][3]<<"), "<< "("<<JQsites[b][4]<<","<<JQsites[b][5]<<"), "<<"("<<JQsites[b][6]<<","<<JQsites[b][7]<<") "; 		
			} 
		} 
		
		
    if (str[0][i]*str[1][i] >=2 ){
			if (k == 0){
					lattice_temp[JHsites[b][0]].flip();
					lattice_temp[JHsites[b][1]].flip();
				}
				else if (k == 1){
					for (int j = 0; j < 2; j++){
						int o = str[j][i];			
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
			//}
		}
		if(str[0][i]*str[1][i] != 0){
				cout<<"   || ms config:  ";
				for (int j = 0; j < Nm; ++j) {
					cout<<lattice_temp[j].S()<<"  ";		
					} 
				cout<<"           "<<endl;
		}					
	}
	for (int i = 0; i < Np; ++i)
				{
					cout << "p: " << i <<" | D/OD: "<<""<<vxoper[pstr[i]]<<" | bnd: ("<<i<<","<<i+Np<<") ";     				
					int o = vxoper[pstr[i]];
					if (o == 2)
					{
						lattice_temp[i].flip();
						lattice_temp[i+Np].flip();
					}
					cout<<"   || ms config:  ";
					for (int j = 0; j < Nm; ++j) {
							cout<<lattice_temp[j].S()<<"  ";		
							} 
					cout<<"           "<<endl;						
				}
}	
void SSEupdates::print_state(SSElattice *lattice)
{
						cout<<" || ms:  " << "  ";
						for (int j = 0; j < Nm; ++j) {
							cout/*<<lattice[i].S()<<"  "*/<<lattice[j].S()<<"  ";		
							}
}					

void SSEupdates::print_op_string(int i)
{
				//cout<<"\noperator string at iteration 'i'= " << w << " \n";
				//for (int i = 0; i < Lc; ++i) {
			 	  //cout << "At slice " << i <<" | I/D/OD is "<<str[0][i]<<","<<str[1][i]<< "  | J/QoQ is "<<str[2][i] <<" | bond loc = "<<str[3][i]-1<<endl;
				//}
				
		int k = str[2][i];
		int b = str[3][i] - 1;

			if (str[0][i]*str[1][i]!=0){
		  if (k==0){
		  cout << "\n";
			cout << "i: " << i <<" | I/D/OD: "<<str[0][i]<< " | J/QoQ: "<<str[2][i] <<" | sites: ("<<JHsites[b][0]<<","<<JHsites[b][1]<<") "; 
			} 
			else if (k == 1){ 
			cout << "\n";
			cout << "i: " << i <<" | I/D/OD: "<<str[0][i]<<","<<str[1][i]<< " | J/QoQ: "<<str[2][i] <<" | sites: ("<<JQsites[b][0]<<","<<JQsites[b][1]<<"), "
			<< "("<<JQsites[b][2]<<","<<JQsites[b][3]<<"), "<< "("<<JQsites[b][4]<<","<<JQsites[b][5]<<"), "<<"("<<JQsites[b][6]<<","<<JQsites[b][7]<<") "; 		
			}
			}
			else{
			cout << "\n";
			cout << "i: " << i <<" | I/D/OD: "<<str[0][i]<< " | J/QoQ: "<<str[2][i] <<" | sites: ("<<0<<","<<0<<") "; 			
			} 
}		
	
	
void SSEupdates::print_proj_op_string(int i)
{
					cout << "\n";
					cout << "p: " << i <<" | D/OD: "<<""<<vxoper[pstr[i]]<<" | bnd: ("<<i<<","<<i+Np<<") ";     						
}

bool SSEupdates::check_states(SSElattice *lattice, int * lattice_prev)
{
			bool check = true;
			for (int j = 0; j < Nm; ++j) {
				if (lattice[j].S() != lattice_prev[j])
				{
					check = false;
					break;
				}	
			}
			return check;
}		

int* SSEupdates::copy_state(SSElattice *lattice)
{
			int* lattice_prev = new int[Nm];
			for (int j = 0; j < Nm; ++j) {
				lattice_prev[j] = lattice[j].S();
			}
			return lattice_prev;
}	
