#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <regex>
#include <iterator>
#include <ctime>
#include <numeric>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "option.h"
#include "pricing_method.h"
#include "Option_price.h"
#include <boost/date_time.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
//#include <bisection.h>

using namespace std;
using namespace boost::gregorian; // compile with linker option as g++ computeBusinessDays.cpp -l boost_date_time

const std::regex comma(",");


long dateDifference( string start_date, string end_date )
{
    date _start_date(from_simple_string(start_date));
    date _end_date(from_simple_string(end_date));

    // counter for weekdays
    int count=0;
    for(day_iterator iter = _start_date; iter!=_end_date; ++iter)
    {
        // increment counter if it's no saturday and no sunday
        if(    iter->day_of_week() !=  boost::date_time::Saturday
            && iter->day_of_week() !=  boost::date_time::Sunday)
            ++count;
    }
    return count;
}

struct Point 
{ 
   string date;
   string exdate;
   string cp_fg;
   double strike_price;
   double  best_bid;
   double  best_offer;

};  

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

double Option_price::getd(){
 return n1;
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
n1 = normalCDF(d1);
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

double Option_price::bisection(double Target,
                 double Low,
                 double High,
                 double Tolerance
                 )
{
    double x=0.5*(Low+High);
    set_sigma(x);
    //cout << x << endl;
    double y=Black_Scholes_Option_Price();
    //cout << y << endl;
 do{
if (y < Target)
    Low = x;
if (y > Target)
    High = x;
x = 0.5*(Low+High);
set_sigma(x);
y = Black_Scholes_Option_Price();
//cout << x << endl;
}while (fabs(y-Target) > Tolerance);
return x; 
}



struct Reader
{
    boost::gregorian::date date1;
    double S;
    double V;
    double ivol;
    double delta;
    double HedE;
    double PNL;
};

double find_Stock(vector<double> st_vec, vector<string> date_vec, string date)
{

    double x;
    int count=1;
   for( int i = 0; i < date_vec.size(); i++)
     {
         if( date == date_vec[i])
         {
             //return rate_vec[i];
             x=st_vec[i];
             break;
         }
        count=count+1;
     }
    if(count > date_vec.size())
    {   
        x=100000.0;
    }

     return x;
}

double find_rate(vector<double> rate_vec, vector<string> date_vec, string date)
{

    double x;
    int count=1;
   for( int i = 0; i < date_vec.size(); i++)
     {
         if( date == date_vec[i])
         {
             //return rate_vec[i];
             x=rate_vec[i];
             break;
         }
        count=count+1;
     }
    if(count > date_vec.size())
    {   
        x=100000.0;
    }

     return x;
}



int main()
{
   vector<string> date;
    vector<double> rate;
    vector<double> Stock_price;
    ifstream infile("./sec_GOOG.csv");
    regex d("^[0-9.-]*$");
    stringstream ss;
    

    string line;
    string w;
    while(!infile.eof())
    {
        getline(infile,line);
        stringstream  ss(line);
        getline(ss ,w,',');
        if (regex_match(w,d))
        {
            date.push_back(w);
            getline(ss,w,',');
            Stock_price.push_back(atof(w.c_str()));
        }
            
    }
   infile.close();
cout << "read interest rate  " << endl;
/*for(auto i = 0; i < rate.size(); i++){
    cout << date[i] << " ";
    cout << rate[i] << endl;
}*/
 
ifstream infile2("./interest.csv");
regex d2("^[0-9.-]*$");
     stringstream ss2;
    

    string line2;
    string w2;
    while(!infile2.eof())
    {
        getline(infile2,line2);
        stringstream  ss2(line2);
        getline(ss2 ,w2,',');
        if (regex_match(w2,d2))
        {
            date.push_back(w2);
            getline(ss2,w2,',');
            rate.push_back(atof(w2.c_str()));
        }
            
    }
   infile2.close();
cout << "read security price  " << endl;
   

 string line1;
    string w1;
regex d1("^[0-9.-]*$");
     stringstream ss1;
Point ass;
vector <Point> sub;  
 ifstream infile1("./op_GOOG.csv");
 while(!infile1.eof())
    {
        getline(infile1,line1);
        stringstream  ss1(line1);
        getline(ss1 ,w1,',');
        if (regex_match(w1,d1))
        {
            //sub.push_back(Point());
            ass.date = w1;
            getline(ss1,w1,',');
            ass.exdate = w1;
            getline(ss1,w1,',');
            ass.cp_fg = w1;
            getline(ss1,w1,',');
            ass.strike_price = stod(w1);
            //sub[0].strike_price = w1;
            getline(ss1,w1,',');
            ass.best_bid = stod(w1);
            getline(ss1,w1,',');
            ass.best_offer = stod(w1);
            //sub[0].best_bid = w1;
             sub.push_back(ass);  
             //cout << ass.date << " " << ass.exdate << " " << ass.cp_fg << " " << ass.strike_price << " " << ass.best_bid << " " << ass.best_offer << endl;                                                                                                                                                                                                                                   //rate.push_back(atof(w1.c_str()));
        }
            
    }
   infile1.close();
   cout << "read option price " << endl;
   

  string mydate;
  string myenddate;
  string my_expdate;
  double K_s;
  
  cout << "please enter a start date" << endl;
  cin >> mydate;
  cout << "please enter an end  date" << endl;
  cin >> myenddate;
  cout << "please enter strike Price" << endl;
  cin  >> K_s;
  cout << "please enter my expire date" << endl;
  cin >> my_expdate;
  long double dd = dateDifference( mydate,my_expdate);
  vector <Reader> filter_data;
  //Reader ass1;
  boost::gregorian::date d2_1 = from_simple_string(mydate);
  long double dd1 = dd;
  
 
   vector <Point> filter_data_is_mine;
for(int i =0;  i < 271890 ; i++)
{
    Point avg;
    if(sub[i].exdate == my_expdate && sub[i].strike_price == K_s && sub[i].cp_fg == "C" ){
        avg.date = sub[i].date;
        avg.exdate = sub[i].exdate;
        avg.cp_fg = sub[i].cp_fg;
        avg.strike_price = sub[i].strike_price;
        avg.best_bid = sub[i].best_bid;
        avg.best_offer = sub[i].best_offer;
        filter_data_is_mine.push_back(avg);
        //cout << avg.date << " " << avg.exdate << " " << avg.cp_fg << " " << avg.strike_price << " " << avg.best_bid << " " << avg.best_offer << endl; 
    }
    //filter_data_is_mine.push_back(avg);
                                                                                                               
}
vector <Point> filter_data_is_mine2;
for(int i = 0; i < filter_data_is_mine.size(); i++){
			if(filter_data_is_mine[i].date <= myenddate && filter_data_is_mine[i].date >= mydate) {
				filter_data_is_mine2.push_back(filter_data_is_mine[i]);			
			}
		}
cout << "filtered data" << endl;
double be=find_Stock(Stock_price, date, mydate);
cout << be << endl;
 double ce = find_rate(rate, date, mydate);
 ce = ce/100.0;
 dd = dd/252;
 cout << ce << endl;
 cout << dd << endl;
  //double sigma =p1.bisection(op,0,1000,0.00001);
	Option_price p0(K_s, be , ce ,dd, 0.01 ,"c");
    double bid = filter_data_is_mine2[0].best_bid;
    double offer = filter_data_is_mine2[0].best_offer;
    double op = 0.5*(bid + offer);
    cout << op << endl;
    double sigma = p0.bisection(op,0.0,20.0,0.00001);
    cout << sigma << endl;
    Option_price p1(K_s, be , ce ,dd, sigma ,"c");
    double opt0 = p1.Black_Scholes_Option_Price();
    cout << opt0 << endl;
   vector <double> tempb;
   vector <double> rate_c;
   for(auto i = 0 ; i < filter_data_is_mine2.size() ; i++)
{  
    cout << "filtered data" << endl;
       mydate = filter_data_is_mine2[i].date;
       d2_1 = from_simple_string(mydate);
      double b = find_Stock(Stock_price, date, mydate);
      double c = find_rate(rate, date, mydate);
       c = c/100.0;
       Option_price p3(K_s, b , c ,dd, sigma ,"c");
      double opt0 = p3.Black_Scholes_Option_Price();
      double bid = filter_data_is_mine2[i].best_bid;
      double offer = filter_data_is_mine2[i].best_offer;
       double op = 0.5*(bid + offer);
      sigma = p3.bisection(op,0.0,20.0,0.01);
      Option_price p4(K_s, b , c ,dd, sigma ,"c");
      double B;
    if(i==0)
       {    
           cout << "filtered data" << endl;
            filter_data.push_back(Reader());
              filter_data[i].date1 = d2_1;
             // cout << filter_data[i].date1 << endl;
              //mydate = to_simple_string(d2_1);
              long  d_e = dateDifference(mydate,my_expdate);
              //cout << d_e << endl;
              //double b = find_Stock(Stock_price, date, mydate);
              //cout << b << endl;
             //double c = find_rate(rate, date, mydate);
             //cout << c << endl;
              //c = c/100.0;
              //dd = dd/252.0;
              be = b;
              //ce = c;
              //d_e = d_e/252.0;
              filter_data[i].S = be;
            
              double del = p3.getd();
              //cout << del << endl;
              double op1 = p3.Black_Scholes_Option_Price();
              double sigma = p3.bisection(op,0.0,20.0,0.01);
              //cout << sigma << endl;
              filter_data[i].V = op1;
              filter_data[i].ivol = sigma;
              filter_data[i].delta = del;
              double pnl = filter_data[0].V;     
              B = op1 - (del*be);
              double H = 0;
              filter_data[i].HedE = H;
              filter_data[i].PNL = pnl;
              tempb.push_back(B);
              rate_c.push_back(c);
         }
    else
         {
             cout << "filtered data" << endl;
              filter_data.push_back(Reader());
              filter_data[i].date1 = d2_1;
              long  d_e = dateDifference(mydate,my_expdate);
              double b = find_Stock(Stock_price, date, mydate);
              double c = find_rate(rate, date, mydate);
              c = c/100.0;
              //dd = dd/252.0;
              be = b;
              ce = c;
              d_e = d_e/252.0;
              filter_data[i].S = be;
              double del = p3.getd();
                //cout << "in here" << endl;double del = p1.getd();
              //cout << del << endl;
              double op1 = p3.Black_Scholes_Option_Price();
              double sigma = p3.bisection(op,0.0,20.0,0.00001);
              //cout << sigma << endl;
              filter_data[i].V = op1;
              filter_data[i].ivol = sigma;
              filter_data[i].delta = del;
              double pnl = filter_data[i].V; 
              double B = (filter_data[i-1].delta*filter_data[i].S) + (tempb[i-1]*exp(rate_c[i-1]*(1/252)))- (filter_data[i].delta*filter_data[i].S);
              double H =  (filter_data[i-1].delta*filter_data[i].S) + (tempb[i-1]*exp(rate_c[i-1]*(1/252))) - filter_data[i].V - filter_data[i-1].V;
              //Delta2[i-1]* sec_price_t+ Bank_balance2[i-1] * exp(rate_t * (RF.dateDifference(filt_vec[i - 1].get_date(), curr_date) / 252.0)) - opt_Price_t
              filter_data[i].HedE = H;
              //cout << filter_data[i].HedE << endl;
              filter_data[i].PNL = pnl;
              //cout << filter_data[i].PNL << endl;
              tempb.push_back(B);
              rate_c.push_back(c);
        }
              
            
           //d2_1 += boost::gregorian::days(1);
           //mydate = to_simple_string(d2_1);

           
}


ofstream fs;
//std::string filename = "exampleOutput.csv";
fs.open("./exampleOutput.csv");
fs << "Date" << "," << "Stock Price" << "," << "Option Price" << "," << "Implied Volatility" << "," << "Delta" << "," << "PNL" << "," << "PNL with Hedge = same as Hedging Error" << std::endl;

for (long i = 0; i < dd1; i++ ){
   fs <<  filter_data[i].date1  << "," << filter_data[i].S  << "," << filter_data[i].V << "," <<  filter_data[i].ivol << "," << filter_data[i].delta << ","  << filter_data[i].PNL <<  ","  << filter_data[i].HedE << std::endl;
   

   //d2_1 += boost::gregorian::days(1);
  }

   



fs.close();



//unit tests
    cout << "testing find stock price function" << endl;
	double unit_test1 = find_Stock(Stock_price, date, "2011-02-24");
	if (unit_test1 == 608.82){
		cout << "unit test passed" << endl;
	}
	else {
		cout << "unit test failed" << endl;
	}

    cout << "testing find rate function" << endl;
	double unit_test2 = find_rate(Stock_price, date, "2011-01-18");
	if (unit_test2 == 0.27 ) {
		cout << "unit test passed" << endl;
	}
	else {
		cout << "unit test failed" << endl;
	}

    cout << "testing bisection function" << endl;
    Option_price p2(100, 100 , 0.05 ,1, 0.01 ,"c");
	double unit_test3 = p2.bisection(10.5,0.05,0.35,0.001);
	if (unit_test3 >= 0.2) {
		cout << "unit test 3 passed" << endl;
	}
	else {
		cout << "unit test 3 failed" << endl;
	}
    
    return 0;
}
