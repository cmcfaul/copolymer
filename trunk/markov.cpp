//$Id$
#include <iostream>
#include <fstream>
#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

int main()
{
    srand(time(0)); // generate seed for random number
    int chains = 10000; // number of polymer chains
    int length = 1000; //length of each polymer chain
    float p_a, p_b; //single-monomer unconditional probabilities
    float p_aa, p_ab, p_ba, p_bb; //diad unconditional probabilities
    float p_aga, p_agb, p_bga, p_bgb; //conditional probabilities
    float r1, r2, x; //input parameters
    int aa, ab, ba, bb; //counters for each diad
    float a_mean, aa_mean, ab_mean, ba_mean, bb_mean;   
    float a_var, aa_var, ab_var, ba_var, bb_var;
    int bad;
    ofstream myfile, myfile2;
    string outfile, statsfile;

    cout << "File name for output ";
    cin >> outfile; statsfile = outfile;
    cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
    cin >> r1 >> r2 >> x;
    
    //calculate the probabilites from the input parameters
    p_a = (x + r1*x*x)/(r1*x*x + 2*x + r2);
    p_b = (x + r2)/(r1*x*x + 2*x + r2);
    p_agb = 1/(1+r2/x);
    p_bga = 1/(1+r1*x);
    p_aga = 1-p_bga;
    p_bgb = 1-p_agb;
    
    //write the header information
    statsfile += "_stats.dat";
    outfile += ".dat";
    myfile.open(outfile.c_str());
    myfile2.open(statsfile.c_str());
    myfile << "Created by $Id$" << endl;
    myfile << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
    myfile2 << "Created by $Id$" << endl;
    myfile2 << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
    myfile << "p_a=" << p_a << " p_b=" << p_b << endl;
    myfile << "p_aga=" << p_aga << " p_agb=" << p_agb << " p_bga=" << p_bga << " p_bgb=" << p_bgb << endl;
    
    //column titles
    myfile << "chain #, a, aa, ab, ba, bb" << endl; 
    //Main body of the program
    a_mean = aa_mean = ab_mean = ba_mean = bb_mean = 0;   
    a_var = aa_var = ab_var = ba_var = bb_var = 0;
    for (int i=0; i<=chains; i++)
    {
	    int monomer = 0;
	    int monomer_a = 0;
	    int monomer_b = 0;
	    aa = ab = ba = bb = 0;
        for (int j=0; j<=length; ++j)
        {
            bad = 0;
	        float k = float(rand())/RAND_MAX;
	        //0 is used to desginate monomer b, 1 means monomer a
	        if (monomer == 0)
               if (k < p_agb)
               {
     		         monomer = 1;
                     monomer_a++;
                     ba++;
                     bad = 1;
               } else
               {
                     monomer = 0;
                     monomer_b++;
                     bb++;
                     bad = 1;
               }
	        else
               if (k < p_aga)
	           {
                     monomer = 1;
                     monomer_a++;
                     aa++;
                     bad = 1;
               } else 
	           {
                     monomer = 0;
   		             monomer_b++;
   		             ab++;
                     bad = 1;
               }
        assert(bad);       
        }
        myfile << i << ", " << monomer_a << ", " << aa << ", " << ab << ", " << ba << ", " << bb << endl;
        a_mean+=monomer_a; aa_mean+=aa; ab_mean+=ab; ba_mean+=ba; bb_mean+=bb;
        a_var+=monomer_a*monomer_a; aa_var+=aa*aa; ab_var+=ab*ab; ba_var+=ba*ba;bb_var+=bb*bb; 
    }
    a_mean /= chains; aa_mean /= chains; ab_mean /= chains; ba_mean /= chains; bb_mean /= chains;
    a_var = a_var/chains - (a_mean*a_mean); aa_var = aa_var/chains - (aa_mean*aa_mean); ab_var = ab_var/chains - (ab_mean*ab_mean);
    ba_var = ba_var/chains - (ba_mean*ba_mean); bb_var = bb_var/chains - (bb_mean*bb_mean);
    myfile2 << "mean, " << a_mean << ", " << aa_mean << ", " << ab_mean << ", " << ba_mean << ", " << bb_mean << endl;
    myfile2 << "std, " << sqrt(a_var) << ", " << sqrt(aa_var) << ", " << sqrt(ab_var) << ", " << sqrt(ba_var) << ", " << sqrt(bb_var) << endl;
    
    myfile.close();
    myfile2.close();
    return 0;
}

