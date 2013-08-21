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

   double ptcut_i = 0.001;
   double ptcut_f = 10;
   int npTcut = 150;
   double dptcut = (log(ptcut_f) - log(ptcut_i))/(npTcut - 1 + 1e-10);

   double *res = new double [2];
   for(int i = 0 ; i < npTcut; i++)
   {
      double ptcut = ptcut_i*exp(i*dptcut);
      double result_tot_eq, result_hard_eq, result_soft_eq;
      double result_tot_vis, result_hard_vis, result_soft_vis;
      //Compton Scattering
      filename = "QGP_2to2_hard";
      test.calculateEmissionrates_hard(filename, ptcut, res);
      result_hard_eq = res[0];
      result_hard_vis = res[1];
   
      //Calculate photon polarization tensor for soft momentum contribution
      filename = "QGP_2to2_soft";
      test.calculateEmissionrates_soft(filename, res);
      result_soft_eq = res[0];
      result_soft_vis = res[1];
      
      result_tot_eq = result_hard_eq + result_soft_eq;
      result_tot_vis = result_hard_vis + result_soft_vis;

      cout << scientific << setprecision(10) << setw(15) 
           << ptcut << "   " << result_hard_eq << "   " 
           << result_soft_eq << "   " << result_tot_eq << "   "
           << result_hard_vis << "   " << result_soft_vis 
           << "   " << result_tot_vis 
           << endl;
   }
   delete [] res;

   sw.toc();
   //cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
