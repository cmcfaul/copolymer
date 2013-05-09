//$Id$
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cmath>

using namespace std;

int main()
{
  int steps = 1000; //number of steps to integrate over
  float h; //time increments
  float x; //beginning ratio of m1/m2
  float r1, r2; //reactivity ratios, input from user
  float m1, m2, dm1, dm2; //actual monomer compositions and differential
  float k11, k12, k21, k22;
  string outfile; 
  ofstream myfile;
  float f1(float m1, float m2, float r1);
  float f2(float m1, float m2, float r2);

  cout << "Name for output file ";
  cin >> outfile; 
  cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
  cin >> r1 >> r2 >> x;

  //convert the input ratio into monomer compositions
  m1 = x/(1.0+x);
  m2 = 1.0/(1.0+x);

  outfile += ".dat";
  myfile.open(outfile.c_str());
  myfile << "Created by $Id$"<< endl;
  myfile << "Time, c_m1, c_m2" << endl;
  for (int i=0; i<=steps; i++)
  {
	  h = float(i)/steps;
	  k11 = h*f1(m1, m2, r1);
	  k12 = h*f2(m1, m2, r2);
	  k21 = h*f1(m1+k11/2, m2+k12/2, r1);
	  k22 = h*f2(m1+k11/2, m2+k12/2, r2);
	  m1 += k21;
	  m2 += k22;
	  assert(m2 >= 0);
	  assert(m1 >= 0);
	  myfile << i << ", " <<  m1 << ", " <<  m2 << endl;

  }
}

float f1(float m1, float m2, float r1)
{
  return -(m1)*(r1*m1 + m2);
}

float f2(float m1, float m2, float r2)
{
  return -(m2)*(r2*m2 + m1);
}

