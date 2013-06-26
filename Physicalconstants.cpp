#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Physicalconstants.h"

Physicalconstants::Physicalconstants()
{
    hbarC = 0.197327053;  //GeV*fm
    e_sq = 1./137.*4.*M_PI;
    N_F = 3.;   // number of flavors
    q_sq = pow(2./3., 2) + pow(-1./3., 2) + pow(-1./3., 2); //u, d, s
    N_c = 3.;   // number of colors
    d_F = 3.;   // dimension of flavor
    C_F = (N_c*N_c - 1)/(2.*N_c);
    alpha_EM = 1./137.;
    g_s = 2.0;

}

Physicalconstants::~Physicalconstants()
{

}


double Physicalconstants::get_g_s_running(double T) 
{
    //Formula from J. Kapusta, "Finite-temperature field theory"
    double g_s = 6.*M_PI/(33. - 2.*N_F)*log(8.*T/0.164);  //T in GeV

    return(g_s);
}
