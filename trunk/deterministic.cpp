//$Id$
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cmath>
//#include "rk4.c" //4th order runge kutta routine from Landau & Paez

using namespace std;

int main()
{
  int steps = 1000; //number of steps to integrate over
  float h; //time increments
  float x; //beginning ratio of m1/m2
  float r[2]; //reactivity ratios, input from user
  float m[2], dm[2]; //actual monomer compositions and differential
  float m0[2];
  float k1[2], k2[2]; 
  string outfile; 
  ofstream myfile;
  float f(float m[2], float r[2], int i);

  cout << "Name for output file ";
  cin >> outfile; 
  cout << "Enter reactivity ratios r1, r2, and input ratio x = [A]/[B]: ";
  cin >> r[0] >> r[1] >> x;

  //convert the input ratio into monomer compositions
  m[0] = m[0] = x/(1.0+x);
  m[1] = m[1] = 1.0/(1.0+x);

  outfile += ".dat";
  myfile.open(outfile.c_str());
  myfile << "Created by $Id$"<< endl;
  myfile << "Time, c_m1, c_m2, f_m1, f_m2, f_total" << endl;
  for (int i=0; i<=steps; i++)
  {
	  h = float(i)/steps;
	  k1[0] = m[0] + 0.5*h*f(m, r, 0);
	  k1[1] = m[1] + 0.5*h*f(m, r, 1);
	  k2[0] = h*f(k1, r, 0);
	  k2[1] = h*f(k1, r, 1);
	  m[0] += k2[0];
	  m[1] += k2[1];
	  assert(m[1] >= 0);
	  assert(m[0] >= 0);
	  myfile << i << ", " <<  m[0] << ", " <<  m[1] << ", ";
	  myfile << (m0[0]-m[0])/m0[0] << ", " << (m0[1]-m[1])/m0[1] << ", " << 1-m[1]-m[0] << endl;

  }
}

float f(float m[2], float r[2], int i)
{
  if (i == 0) return -(m[0])*(r[0]*m[0] + m[1]);
  if (i == 1) return -(m[1])*(r[1]*m[1] + m[0]);
}
