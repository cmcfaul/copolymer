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
  int state; //state of the chain-end
  float a_mean, aa_mean, ab_mean, ba_mean, bb_mean;   
  float a_var, aa_var, ab_var, ba_var, bb_var;
  ofstream myfile, myfile2;
  string outfile, statsfile;

  cout << "Base name for output ";
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
  p_aa = p_aga*p_a;
  p_bb = p_bgb*p_b;
  p_ab = p_bga*p_a;
  p_ba = p_agb*p_b;

  //write the header information
  outfile += ".dat";
  statsfile += "_stats.dat";
  myfile.open(outfile.c_str());
  myfile2.open(statsfile.c_str());
  myfile << "Created by $Id$" << endl;
  myfile << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
  myfile << "p_a=" << p_a << " p_b=" << p_b << endl;
  myfile << "p_aga=" << p_aga << " p_agb=" << p_agb << " p_bga=" << p_bga << " p_bgb=" << p_bgb << endl;
  myfile << "p_aa=" << p_aa << " p_bb=" << p_bb << " p_ab=" << p_ab << " p_ba=" << p_ba << endl;
  myfile2 << "Created by $Id$" << endl;
  myfile2 << "r1=" << r1 << " r2=" << r2 << " [A]/[B]=" << x << endl;
  //column titles
  myfile << "chain #, a, aa, ab, ba, bb, a_left, b_left" << endl;

  //Main body of the program
  a_mean = aa_mean = ab_mean = ba_mean = bb_mean = 0;   
  a_var = aa_var = ab_var = ba_var = bb_var = 0;
  int a_left = int(floor(chains*length*x/(x+1))); //=total monomers * fraction of a monomers
  int b_left = int(floor(chains*length*1/(x+1))); //=total monomers * fraction of b monomers
  for (int i=1; i<=chains; i++)
  {
	  //start chain using initial probabilities
	  float k = float(rand())/(RAND_MAX + 1.0);
	  int a = 0;
	  int b = 0;
	  aa = ab = ba = bb = 0;
	  if (k < p_a)
	  {
		  state = 1;
		  a++;
		  a_left--;
	  }
	  else 
	  {
		  state = 0;
		  b++;
		  b_left--;
	  }
	  for (int j=1; j<length; ++j)
	  {
		  k = float(rand())/(RAND_MAX + 1.0);
		  if (k < p_a)
		  {
			  (state) ? (aa++) : (ba++);
			  state = 1;
			  a++;
			  a_left--;
		  }
		  else
		  {
			  (state) ? (ab++) : (bb++);
			  state = 0;
			  b++;
			  b_left--;
		  }
		  //update probabilities
		  if (a_left == 0 && b_left != 0)
		  {
			  p_a = 0;
			  p_b = 1;
		  }
		  else if (a_left != 0 && b_left == 0)
		  {
			  p_a = 1;
			  p_b = 0;
		  }
		  else if (a_left == 0 && b_left == 0)
		  {
			  //we should be on the last loop
			  assert (i == chains);   
			  assert (j == length-1);
		  }
		  else if (state)//a_left, b_left > 0, state = 1 
		  {
			  x = a_left/b_left;
			  p_b = 1/(1+r1*x);
			  p_a = 1-p_b;
		  }
		  else //a_left, b_left > 0, state = 0
		  {
			  x = a_left/b_left;
			  p_a = 1/(1+r2/x);
			  p_b = 1-p_a;
		  }
	  } //next monomer
	  myfile << i << ", " << a << ", " << aa << ", " << ab << ", " << ba << ", " << bb << ", " << a_left << ", " << b_left << endl;
	  a_mean+=a; aa_mean+=aa; ab_mean+=ab; ba_mean+=ba; bb_mean+=bb;
	  a_var+=a*a; aa_var+=aa*aa; ab_var+=ab*ab; ba_var+=ba*ba; bb_var+=bb*bb; 
  } //next chain
  assert(a_left == 0);
  assert(b_left == 0);

  //compute statistics for the reaction
  a_mean /= chains; aa_mean /= chains; ab_mean /= chains; ba_mean /= chains; bb_mean /= chains;
  a_var = a_var/chains - (a_mean*a_mean); aa_var = aa_var/chains - (aa_mean*aa_mean); ab_var = ab_var/chains - (ab_mean*ab_mean);
  ba_var = ba_var/chains - (ba_mean*ba_mean); bb_var = bb_var/chains - (bb_mean*bb_mean);
  myfile2 << "statistic, a, aa, ab, ba, bb" << endl;
  myfile2 << "mean, " << a_mean << ", " << aa_mean << ", " << ab_mean << ", " << ba_mean << ", " << bb_mean << endl;
  myfile2 << "std, " << sqrt(a_var) << ", " << sqrt(aa_var) << ", " << sqrt(ab_var) << ", " << sqrt(ba_var) << ", " << sqrt(bb_var) << endl;

  myfile.close();
  myfile2.close();
  return 0;
}
