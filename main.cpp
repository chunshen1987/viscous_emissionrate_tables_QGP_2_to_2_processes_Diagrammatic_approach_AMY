#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Arsenal.h"
#include "Stopwatch.h"
#include "QGP_2to2_Scattering.h"
#include "ParameterReader.h"

using namespace std;

int main(int argc, char** argv)
{
   Stopwatch sw; 
   sw.tic();
   
   string filename;
   
   ParameterReader* paraRdr = new ParameterReader();
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);

   QGP_2to2_Scattering test(paraRdr);

   double ptcut_i = 0.005;
   double ptcut_f = 100;
   int npTcut = 30;
   double dptcut = (log(ptcut_f) - log(ptcut_i))/(npTcut - 1);

   for(int i = 0 ; i < npTcut; i++)
   {
      double ptcut = ptcut_i*exp(i*dptcut);
      double result_tot, result_hard, result_soft;
      //Compton Scattering
      filename = "QGP_2to2_hard";
      result_hard = test.calculateEmissionrates_hard(filename, ptcut);
   
      //Calculate photon polarization tensor for soft momentum contribution
      filename = "QGP_2to2_soft";
      result_soft = test.calculateEmissionrates_soft(filename);
      
      result_tot = result_hard + result_soft;
      cout << scientific << setprecision(10) << setw(15) 
           << ptcut << "   " << result_hard << "   " 
           << result_soft << "   " << result_tot 
           << endl;
   }

   sw.toc();
   //cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
