#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <ctime>
#include <numeric>
using namespace std;

double average(vector<double> v)
{
// code for calculating average of members of v
// and returning the average
double sum_of_v;
int no_of_v;
double av;
sum_of_v = accumulate(v.begin(),v.end(),0.0 );
no_of_v = v.size();
av = (sum_of_v/no_of_v);
return av;

}

double find_rate(vector<double> rate_vec, vector<string> date_vec, double date)
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
    stringstream ss;
    ifstream infile("./hw1_H.15_Baa_Data.csv");
    regex d("^[0-9.-]*$");
   // regex r("^[0-9]+(\.[0-9]{1,2})?$");
   
// code for loading rate and date vectors from the file hw1_H.15_Baa_Data.csv
// the headers should be handled properly. do not delete them manually
  /*if (infile.fail()) {
        cerr << "Error opening a file" << endl;
        infile.close();
        exit(1);
    }*/

    string line;
    string w;
    while(!infile.eof())
    {
        getline(infile,line);
        stringstream  ss(line);
        getline(ss ,w,',');
        if (regex_match(w,d))
        {
            year.push_back(w);
            getline(ss,w,',');
            rate.push_back(atof(w.c_str()));
        }
            
    }
   infile.close();
   

// code for prompting user for a date and returning the rate
// and the difference between the rate for that date and the
// average rate
//
string s;
 cout << "please enter a date" << endl;
 while(getline(cin,s)){
 
 if(regex_match(s,d))
 {
    double a = average(rate);
    double b = find_rate(rate, date, s);
    double c = b - a;
    if (b == 100000.0)
    {
        cout << "date not found" << endl;
    }
   else 
   {
       cout << b << " " << c << endl;
   } 
    cout <<  "Enter one more date " << endl;

 }

else 
{
  cout << "Please enter a valid value" << endl;
}
if(cin.eof())

{
    cout << " the end ";
}



 }
// This code should allow the user to continue to input dates
// until the user inputs the EOF (End-of-file), namely control-d (in linux/mac) or control-z (in windows)
// This should not crash if a bad date is entered.

   return 0; // program end
}
