//$Id$
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

int main()
{
  float x; //beginning ratio of m1/m2
  float r1, r2; //reactivity ratios, input from user
  float m1, m2, dm1, dm2; //actual monomer compositions and differential
  float m10, m20; //initial compositions
  int steps = 100; //number of steps to integrate over
  string outfile; 
  ofstream myfile;

  cout << "Name for output file ";
  cin >> outfile; 
  cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
  cin >> r1 >> r2 >> x;

  //convert the input ratio into monomer compositions
  m1 = m10 = x/(1.0+x);
  m2 = m20 = 1.0/(1.0+x);

  myfile.open(outfile.c_str());
  myfile << "Time, c_m1, c_m2" << endl;
  for (int i=0; i<=steps; i++)
  {
	  dm2 = -m2*float(i)/steps;
	  m2 += dm2; 
	  dm1 = dm2*x*(r1*m1 + m2)/(m1 + r2*m2);
	  m1 += dm1;
	  assert(m2 >= 0);
	  assert(m1 >= 0);
	  myfile << i << ", " <<  m1 << ", " <<  m2 << endl;

  }
}

