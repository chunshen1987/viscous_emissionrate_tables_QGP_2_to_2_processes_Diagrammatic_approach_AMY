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
   test.calculateEmissionrates_AMY();

   ////Compton Scattering
   //filename = "QGP_2to2_hard";
   //test.calculateEmissionrates_hard(filename);
   //
   ////Calculate photon polarization tensor for soft momentum contribution
   //filename = "QGP_2to2_soft";
   //test.calculateEmissionrates_soft(filename);

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
