//$Id$
#include <iostream>
#include <fstream>
#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()

using namespace std;

int main()
{
    srand(time(0)); // generate seed for random number
    int chains = 10000; // number of polymer chains
    int length = 1000; //length of each polymer chain
    float p_a, p_b; //single-monomer unconditional probabilities
    float p_aa, p_ab, p_ba, p_bb; //diad unconditional probabilities
    float p_aga, p_agb, p_bga, p_bgb; //conditional probabilities
    float r1, r2, x;
    ofstream myfile;
    char* outfile;

    cout << "File name for output ";
    cin >> outfile;
    cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
    cin >> r1 >> r2 >> x;
    
    //calculate the probabilites from the input parameters
    p_a = (x + r1*x*x)/(r1*x*x + 2*x + r2);
    p_b = (x + r2)/(r1*x*x + 2*x + r2);
    p_agb = 1/(1+r2);
    p_bga = 1/(1+r1);
    p_aga = 1-p_bga;
    p_bgb = 1-p_agb;
    
    //write the header information
    myfile.open(outfile);
    myfile << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
    myfile << "p_a=" << p_a << " p_b=" << p_b << endl;
    myfile << "p_aga=" << p_aga << " p_agb=" << p_agb << " p_bga=" << p_bga << " p_bgb=" << p_bgb << endl;
    myfile << "[A]/[B] = " << x << endl;
    
    for (int i=0; i<=chains; i++)
    {
	int monomer = 0;
	int monomer_a = 0;
        for (int j=0; j<=length; ++j)
	{
//	    myfile << monomer << ", ";
	    float k = float(rand())/RAND_MAX;
	    if (monomer == 0 && k < p_agb)
	    {
     		monomer = 1;
            monomer_a++;;
	    }
	    else if (monomer == 1 && k < p_aga)
	    {
		    monomer = 1;
            monomer_a++;;
	    }
	    else
	    {
   		    monomer = 0;
        }

//	myfile << monomer << endl;
	}
	myfile << monomer_a << endl;
    }
    myfile.close();
	
    return 0;
}

