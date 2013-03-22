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
    float p_A, p_B;
    float r1, r2, x;
    ofstream myfile;
    char* outfile;

    cout << "File name for output ";
    cin >> outfile;
    cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
    cin >> r1 >> r2 >> x;
    p_A = (x + r1*x*x)/(r1*x*x + 2*x + r2);
    p_B = (x + r2)/(r1*x*x + 2*x + r2);
    myfile.open(outfile);
    myfile << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
    myfile << "p_A=" << p_A << " p_B=" << p_B << endl;
    myfile << "[A]/[B] = " << x << endl;
    for (int i=0; i<=chains; i++)
    {
	int monomer = 0;
        for (int j=0; j<=length; ++j)
	{
	   float k = float(rand())/RAND_MAX;
	    if (k < p_A)
	    {
	    monomer++;
	    }
	}
	myfile << monomer << endl;
    }
    myfile.close();
	
    return 0;
}

