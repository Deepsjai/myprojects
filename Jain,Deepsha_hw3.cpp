#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctime>
#include <numeric>
#include <string>
#include <iomanip>
//#include <regex>
#include "option.h"
#include "pricing_method.h"
#include "Option_price.h"
using namespace std;

void option::init()									//init function
{
	StikePrice = 0.0;
	CPU = 0.0;
	r = 0.0;
    t =0.0;
    sigma = 0.0;
}

option::option()								//Default Constructor
{
   init();
}

option::option(double sp, double cpu, double b, double time, double sig)		//Another Constructor
{
	StikePrice = sp;
	CPU = cpu;
	r = b;
    t = time;
    sigma = sig;
}

option::~option()								//Destructor
{
}


Option_price::Option_price(double sp, double cpu, double b, double time, double sig, string f):option(sp, cpu, b, time,sig){
    flag = f;                       // constructor for the flag
}
Option_price::~Option_price()								//Destructor
{
}

double Option_price::Binomial_Option_Price()
{
  double r = get_r();
  double K = get_stikeprice();
  double sigma= get_sigma();
  double S = get_CPU();
  int n_p = 200;
  
  // no. of time periods
  double time = get_t(); //
  double r_inv = exp(-r*(time/n_p));
  double u = exp(sigma* sqrt(time/n_p));
  double d = 1/u;
  double p_up = (exp(r*(time/n_p))-d)/(u-d);
  double p_d = 1.0 - p_up;
  double u_u = u*u;

 std::vector<double> prices(n_p+1); 
  prices[0] = S*pow(d, n_p);// finding the lowest node element of the tree
  for (int i=1; i<=n_p; ++i) prices[i] = u_u*prices[i-1];
  std::vector<double> c_values(n_p+1); // value of corresponding call
  for (int i=0; i<=n_p; ++i) { 
      if(flag == "c"|| flag == "C")c_values[i] = fmax(0.0, (prices[i]-K));// flag for call option
      if(flag == "p"|| flag == "P")c_values[i] = fmax(0.0, (K-prices[i]));//flag for put
      }; //  payoffs at maturity 
for (int step=n_p-1; step>=0; --step) {
  for (int i=0; i<=step; ++i) {
  c_values[i] = (p_up*c_values[i+1]+ p_d*c_values[i])*r_inv;
}; 
};
return c_values[0]; 
}

double normalCDF(double x) 
{
    return erfc(-x / sqrt(2.0))/2.0;
}


double Option_price::Black_Scholes_Option_Price(){
double K = get_stikeprice();
double S = get_CPU();
double r = get_r();
double time = get_t();
double sigma = get_sigma();
double timesqrt = sqrt(time);
double d1 = (log(S/K)+r*time)/(sigma*timesqrt)+0.5*sigma*timesqrt; 
double d2 = d1-(sigma*timesqrt);
double n1 = normalCDF(d1);
double n2 = normalCDF(d2);
double c_price;
if (flag == "c"|| flag == "C"){
c_price = S*n1 - K*exp(-r*time)*n2; 
}
if (flag == "p"|| flag == "P"){
c_price = K*exp(-r*time)*(normalCDF(-d2)) - S*(normalCDF(-d1)); 
}
return c_price;
}



int main()
{

//code for initialising option : Q1 

double x,y,z, v,w;
string c_p;
std::cout << " please enter values for StrikePrice:" << endl;
std::cin >> x;
std::cout << " please enter values for Current price of underlying:" << endl;
std::cin >> y;
std::cout << " please enter values for risk free rate" << endl;
std::cin >> z;
std::cout << " please enter values for time" << endl;
std::cin >> v;
std::cout << " please enter values for volatility" << endl;
std::cin >> w;
std::cout << " please enter values for flag" << endl;
std::cin >> c_p;
option o1(x,y,z,v,w);
double o_x = o1.get_CPU();
std::cout << o_x << endl;
//double s = 125, r = 0.025, u = 1.05, d =1/u, p = 0.745;
//int n_p = 2;
//string flag = "c";
Option_price p1(x,y,z,v,w,c_p);
double C_bin = p1.Binomial_Option_Price();
double C_black =p1.Black_Scholes_Option_Price();
std::cout << "Binomial option price:" << C_bin << " " << endl;
std::cout << "Black_Scholes option price:" << C_black << " " << endl;
return 0;

}
