#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "QGP_2to2_Scattering.h"
#include "Physicalconstants.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"
#include "Arsenal.h"

using namespace std;

QGP_2to2_Scattering::QGP_2to2_Scattering(ParameterReader* paraRdr_in)
{
   paraRdr = paraRdr_in;
   
   n_Eq = paraRdr->getVal("n_Eq");
   double Eq_i = paraRdr->getVal("Eq_min");
   double dEq = paraRdr->getVal("dEq");
   Eq_tb = new double [n_Eq];
   for(int i=0; i<n_Eq; i++)
      Eq_tb[i] = Eq_i + i*dEq;

   n_Temp = paraRdr->getVal("n_Temp");
   double T_i = paraRdr->getVal("T_min");
   double dT = paraRdr->getVal("dT");
   T_tb = new double [n_Temp];
   for(int i=0; i<n_Temp; i++)
      T_tb[i] = T_i + i*dT;

   equilibrium_results = new double* [n_Eq];
   viscous_results = new double*[n_Eq];
   for(int i=0; i<n_Eq; i++)
   {
      equilibrium_results[i] = new double [n_Temp];
      viscous_results[i] = new double [n_Temp];
   }

   n_ktilde = paraRdr->getVal("n_ktilde");
   double ktilde_i = paraRdr->getVal("ktilde_i");
   double ktilde_f = paraRdr->getVal("ktilde_max");
   double dktilde = (ktilde_f - ktilde_i)/(n_ktilde - 1 + 1e-10);
   ktilde_pt = new double [n_ktilde];
   for(int i=0; i<n_ktilde; i++)
      ktilde_pt[i] = ktilde_i + i*dktilde;
   equilibriumTilde_results = new double [n_ktilde];
   viscousTilde_results = new double [n_ktilde];
   
   //initialize the Gaussian quadrature lattices
   n_qtilde = paraRdr->getVal("n_qtilde");
   qtilde_pt = new double [n_qtilde];
   qtilde_weight = new double [n_qtilde];
   qtilde_pt_standard = new double [n_qtilde];
   qtilde_weight_standard = new double [n_qtilde];

   n_omega = paraRdr->getVal("n_omega");
   omega_pt = new double [n_omega];
   omega_weight = new double [n_omega];
   omega_pt_standard = new double [n_omega];
   omega_weight_standard = new double [n_omega];

   n_pprime = paraRdr->getVal("n_pprime");
   pprime_max = paraRdr->getVal("pprime_max");
   pprime_pt = new double [n_pprime];
   pprime_weight = new double [n_pprime];
   pprime_pt_standard = new double [n_pprime];
   pprime_weight_standard = new double [n_pprime];

   n_pSoft = paraRdr->getVal("n_p");
   pSoft = new double [n_pSoft];
   pSoft_weight = new double [n_pSoft];

   n_theta = paraRdr->getVal("n_theta");
   costhetaSoft = new double [n_theta];
   costhetaSoft_weight = new double [n_theta];

   deltaf_alpha = paraRdr->getVal("deltaf_alpha");

   set_gausspoints();
   
}

QGP_2to2_Scattering::~QGP_2to2_Scattering()
{
   delete[] Eq_tb;
   delete[] T_tb;
   for(int i=0; i<n_Eq; i++)
   {
      delete[] equilibrium_results[i];
      delete[] viscous_results[i];
   }
   delete[] equilibrium_results;
   delete[] viscous_results;
   delete[] equilibriumTilde_results;
   delete[] viscousTilde_results;

   delete[] ktilde_pt;

   delete[] qtilde_pt;
   delete[] qtilde_weight;
   delete[] qtilde_pt_standard;
   delete[] qtilde_weight_standard;
   delete[] omega_pt;
   delete[] omega_weight;
   delete[] omega_pt_standard;
   delete[] omega_weight_standard;
   delete[] pprime_pt;
   delete[] pprime_weight;
   delete[] pprime_pt_standard;
   delete[] pprime_weight_standard;

   delete [] pSoft;
   delete [] pSoft_weight;
   delete [] costhetaSoft;
   delete [] costhetaSoft_weight;
}

void QGP_2to2_Scattering::set_gausspoints()
{
   gauss_quadrature_standard(n_qtilde, 1, 0.0, 0.0, 0.0, 1.0, qtilde_pt_standard, qtilde_weight_standard);

   gauss_quadrature_standard(n_omega, 1, 0.0, 0.0, 0.0, 1.0, omega_pt_standard, omega_weight_standard);
   
   gauss_quadrature_standard(n_pprime, 1, 0.0, 0.0, 0.0, 1.0, pprime_pt_standard, pprime_weight_standard);
  
   gauss_quadrature(n_theta, 1, 0.0, 0.0, -1., 1., costhetaSoft, costhetaSoft_weight);

}

void QGP_2to2_Scattering::buildupEmissionrate2DTable()
{
   double *ktildeT = new double [n_Eq];
   double *rawResult_eq = new double [n_Eq];
   double *rawResult_vis = new double [n_Eq];
   double *log_eq = new double [n_ktilde];
   double *scaled_vis = new double [n_ktilde];
   for(int i = 0; i < n_ktilde; i++)
   {
      log_eq[i] = log(equilibriumTilde_results[i]);
      scaled_vis[i] = viscousTilde_results[i]/equilibriumTilde_results[i];
   }
   for(int j = 0; j < n_Temp; j++)
   {
      double T = T_tb[j];
      for(int i = 0; i < n_Eq; i++)
         ktildeT[i] = Eq_tb[i]/T;
      interpolation1D_linear(ktilde_pt, log_eq, ktildeT, rawResult_eq, n_ktilde, n_Eq);
      interpolation1D_linear(ktilde_pt, scaled_vis, ktildeT, rawResult_vis, n_ktilde, n_Eq);
      for(int i = 0; i < n_Eq; i++)
      {
         double temp = exp(rawResult_eq[i]);
         equilibrium_results[i][j] = T*T*temp;
         viscous_results[i][j] = rawResult_vis[i]*temp;
      }
   }
   delete [] ktildeT;
   delete [] rawResult_eq;
   delete [] rawResult_vis;
   delete [] log_eq;
   delete [] scaled_vis;
}

void QGP_2to2_Scattering::output_emissionrateTable()
{
   // output emission rate tables
   ostringstream output_file_eqrate;
   ostringstream output_file_viscous;
   output_file_eqrate << "rate_" << filename << "_eqrate.dat";
   output_file_viscous << "rate_" << filename << "_viscous.dat";
   ofstream of_eqrate(output_file_eqrate.str().c_str());
   ofstream of_viscous(output_file_viscous.str().c_str());
   for(int j=0; j<n_Temp; j++)
   {
      for(int i=0; i<n_Eq; i++)
      {
         of_eqrate << scientific << setw(20) << setprecision(8)
                   << equilibrium_results[i][j] << "   ";
         of_viscous << scientific << setw(20) << setprecision(8)
                   << viscous_results[i][j] << "   ";
      }
      of_eqrate << endl;
      of_viscous << endl;
   }

   of_eqrate.close();
   of_viscous.close();
}

void QGP_2to2_Scattering::calculateEmissionrates_I1(int channel_in, string filename_in)
{
   double hbarC = Phycons.get_hbarC();
   double e_sq = Phycons.get_e_sq();
   double q_sq = Phycons.get_q_sq();
   double d_F = Phycons.get_d_F();
   double C_F = Phycons.get_C_F();
   double g_s = Phycons.get_g_s_const();
   channel = channel_in;
   filename = filename_in;

   double *results = new double [2];
   for(int i=0; i<n_ktilde; i++)
   {
       double ktilde = ktilde_pt[i];

       double prefactor = 1./(pow(2.0*M_PI, 6)*ktilde)*e_sq*q_sq*d_F*C_F*g_s*g_s;

       scale_gausspoints_omega(ktilde);

       double equilibrium_result_omega = 0.0;
       double viscous_result_omega = 0.0;

       for(int k=0; k<n_omega; k++)
       {
          Integrate_I1_qtilde(ktilde, omega_pt[k], results);

          equilibrium_result_omega += results[0]*omega_weight[k];
          viscous_result_omega += results[1]*omega_weight[k];
       }
       equilibriumTilde_results[i] = equilibrium_result_omega*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
       viscousTilde_results[i] = viscous_result_omega*prefactor/(ktilde*ktilde)/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
       
       cout << equilibriumTilde_results[i]*pow(hbarC,4) << endl;
       //cout << viscousTilde_results[i] << endl;

   }
   
   delete [] results;
}

void QGP_2to2_Scattering::scale_gausspoints_omega(double ktilde)
{
   double omega_min = ktilde;
   double omega_max = paraRdr->getVal("omega_max");
   if( ktilde*200 > omega_max) omega_max = ktilde*200;
  
   for(int i=0; i < n_omega; i++)
   {
      omega_pt[i] = omega_pt_standard[i];
      omega_weight[i] = omega_weight_standard[i];
   }
   scale_gausspoints(n_omega, 1, 0.0, 0.0, omega_min, omega_max, omega_pt, omega_weight);
  
}

void QGP_2to2_Scattering::Integrate_I1_qtilde(double ktilde, double omega, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;

   double qtilde_min = fabs(2.*ktilde - omega);
   double qtilde_max = omega;
   for(int i = 0; i < n_qtilde; i++)
   {
      qtilde_pt[i] = qtilde_pt_standard[i];
      qtilde_weight[i] = qtilde_weight_standard[i];
   }
   scale_gausspoints(n_qtilde, 1, 0.0, 0.0, qtilde_min, qtilde_max, qtilde_pt, qtilde_weight);
   for(int i = 0; i < n_qtilde; i++)
   {
      Integrate_I1_pprime(ktilde, omega, qtilde_pt[i], results);
      equilibrium_result += results[0]*qtilde_weight[i];
      viscous_result += results[1]*qtilde_weight[i];
   }

   results[0] = equilibrium_result;
   results[1] = viscous_result;
}

void QGP_2to2_Scattering::Integrate_I1_pprime(double ktilde, double omega, double qtilde, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;

   double pprime_min = (omega - qtilde)/2.;
   double pprime_max_I1 = (omega + qtilde)/2.;
   for(int i = 0; i < n_pprime; i++)
   {
      pprime_pt[i] = pprime_pt_standard[i];
      pprime_weight[i] = pprime_weight_standard[i];
   }
   scale_gausspoints(n_pprime, 1, 0.0, 0.0, pprime_min, pprime_max_I1, pprime_pt, pprime_weight);

   double manderstan_s = omega*omega - qtilde*qtilde;
   double cos_theta_kq = (qtilde*qtilde - omega*omega + 2.*omega*ktilde)/(2.*qtilde*ktilde);
   double sin_theta_kq = sqrt(1. - cos_theta_kq*cos_theta_kq);

   // - t/s term
   for(int i = 0; i < n_pprime; i++)
   {
      double pprime = pprime_pt[i];
      double cos_theta_pprimeq = (qtilde*qtilde - omega*omega + 2.*omega*pprime)/(2.*pprime*qtilde);
      double sin_theta_pprimeq = sqrt(1. - cos_theta_pprimeq*cos_theta_pprimeq);
      double f0_E1 = Bose_distribution(omega-pprime);
      double f0_E2 = Fermi_distribution(pprime);
      double f0_E3 = Fermi_distribution(omega-ktilde);
      double eq_integrand = 2.*ktilde/manderstan_s*((omega - pprime) - cos_theta_kq*(qtilde - pprime*cos_theta_pprimeq))*f0_E1*f0_E2*(1. - f0_E3);
      double vis_integrand = eq_integrand*((1. + f0_E1)*deltaf_chi(omega-pprime)*(-0.5 + 1.5*Power(1./(omega - pprime), 2)*(qtilde*qtilde*cos_theta_kq*cos_theta_kq + 0.5*pprime*pprime*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq + pprime*pprime*cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 2.*qtilde*pprime*cos_theta_pprimeq*cos_theta_kq*cos_theta_kq) )
                                          +(1. - f0_E2)*deltaf_chi(pprime)*(-0.5 + 1.5*(cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 0.5*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq))
                                          -f0_E3*deltaf_chi(omega-ktilde)*(-0.5 + 1.5*Power((qtilde*cos_theta_kq - ktilde)/(omega - ktilde), 2))
                                          )
                            - pprime*ktilde/manderstan_s*sin_theta_kq*sin_theta_pprimeq*f0_E1*f0_E2*(1. - f0_E3)
                             *( (1. - f0_E2)*deltaf_chi(pprime)*(3.*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq*cos_theta_pprimeq)
                               +(1. + f0_E1)*deltaf_chi(omega - pprime)*(3.*Power(1./(omega - pprime), 2)*(- pprime*qtilde*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq - pprime*pprime*sin_theta_pprimeq*cos_theta_pprimeq*sin_theta_kq*cos_theta_kq))
                              );
      equilibrium_result += eq_integrand*pprime_weight[i];
      viscous_result += vis_integrand*pprime_weight[i];
   }

   results[0] = equilibrium_result;
   results[1] = viscous_result;
}

void QGP_2to2_Scattering::calculateEmissionrates_I2(int channel_in, string filename_in)
{
   double hbarC = Phycons.get_hbarC();
   double e_sq = Phycons.get_e_sq();
   double q_sq = Phycons.get_q_sq();
   double d_F = Phycons.get_d_F();
   double C_F = Phycons.get_C_F();
   double g_s = Phycons.get_g_s_const();
   channel = channel_in;
   filename = filename_in;

   double *results = new double [2];
   for(int i=0; i<n_ktilde; i++)
   {
       double ktilde = ktilde_pt[i];

       double prefactor = 1./(pow(2.0*M_PI, 6)*ktilde)*e_sq*q_sq*d_F*C_F*g_s*g_s;

       scale_gausspoints_qtilde(ktilde);

       double equilibrium_result_qtilde = 0.0;
       double viscous_result_qtilde = 0.0;

       for(int k=0; k<n_qtilde; k++)
       {
          Integrate_I2_omega(ktilde, qtilde_pt[k], results);

          equilibrium_result_qtilde += results[0]*qtilde_weight[k];
          viscous_result_qtilde += results[1]*qtilde_weight[k];
       }
       equilibriumTilde_results[i] = equilibrium_result_qtilde*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
       viscousTilde_results[i] = viscous_result_qtilde*prefactor/(ktilde*ktilde)/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
       
       cout << equilibriumTilde_results[i]*pow(hbarC,4) << endl;
       //cout << viscousTilde_results[i] << endl;
   }
   
   delete [] results;
}


void QGP_2to2_Scattering::scale_gausspoints_qtilde(double ktilde)
{
   double g_s = Phycons.get_g_s_const();
   qtilde_cutoff = sqrt(g_s);
   double qtilde_min = qtilde_cutoff;
   double qtilde_max = paraRdr->getVal("qtilde_max");
   if( ktilde*200 > qtilde_max) qtilde_max = ktilde*200;
  
   for(int i=0; i < n_qtilde; i++)
   {
      qtilde_pt[i] = qtilde_pt_standard[i];
      qtilde_weight[i] = qtilde_weight_standard[i];
   }
   scale_gausspoints(n_qtilde, 1, 0.0, 0.0, qtilde_min, qtilde_max, qtilde_pt, qtilde_weight);
  
}

void QGP_2to2_Scattering::Integrate_I2_omega(double ktilde, double qtilde, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;

   double omega_min = qtilde - 2*ktilde;
   if(omega_min < - qtilde)
      omega_min = - qtilde;
   double omega_max = qtilde;
   for(int i = 0; i < n_omega; i++)
   {
      omega_pt[i] = omega_pt_standard[i];
      omega_weight[i] = omega_weight_standard[i];
   }
   scale_gausspoints(n_omega, 1, 0.0, 0.0, omega_min, omega_max, omega_pt, omega_weight);
   for(int i = 0; i < n_omega; i++)
   {
      Integrate_I2_pprime(ktilde, qtilde, omega_pt[i], results);
      equilibrium_result += results[0]*omega_weight[i];
      viscous_result += results[1]*omega_weight[i];
   }

   results[0] = equilibrium_result;
   results[1] = viscous_result;
}

void QGP_2to2_Scattering::Integrate_I2_pprime(double ktilde, double qtilde, double omega, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;

   double pprime_min = (qtilde - omega)/2.;
   for(int i = 0; i < n_pprime; i++)
   {
      pprime_pt[i] = pprime_pt_standard[i];
      pprime_weight[i] = pprime_weight_standard[i];
   }
   scale_gausspoints(n_pprime, 1, 0.0, 0.0, pprime_min, pprime_max, pprime_pt, pprime_weight);

   double manderstan_t = omega*omega - qtilde*qtilde;
   double cos_theta_kq = (omega*omega - qtilde*qtilde + 2.*omega*ktilde)/(2.*qtilde*ktilde);
   double sin_theta_kq = sqrt(1. - cos_theta_kq*cos_theta_kq);

   double term_st, term_ut, term_st_vis, term_ut_vis;
   // - s/t term
   for(int i = 0; i < n_pprime; i++)
   {
      double pprime = pprime_pt[i];
      double cos_theta_pprimeq = (omega*omega - qtilde*qtilde + 2.*omega*pprime)/(2.*pprime*qtilde);
      double sin_theta_pprimeq = sqrt(1. - cos_theta_pprimeq*cos_theta_pprimeq);
      double f0_E1 = Bose_distribution(omega+ktilde);
      double f0_E2 = Fermi_distribution(pprime);
      double f0_E3 = Fermi_distribution(omega+pprime);
      double eq_integrand = - (-1. - (- 2.*pprime*ktilde)/manderstan_t*(1. - cos_theta_pprimeq*cos_theta_kq))*f0_E1*f0_E2*(1. - f0_E3);
      double vis_integrand = eq_integrand*((1. - f0_E1)*deltaf_chi(omega+ktilde)*(-0.5 + 1.5*Power((qtilde*cos_theta_kq + ktilde)/(omega + ktilde), 2))
                                          +(1. - f0_E2)*deltaf_chi(pprime)*(-0.5 + 1.5*(cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 0.5*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq))
                                          +f0_E3*deltaf_chi(pprime+omega)*(-0.5 + 1.5*Power(1./(pprime+omega), 2)*(0.5*pprime*pprime*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq + qtilde*qtilde*cos_theta_kq*cos_theta_kq + pprime*pprime*cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 2.*qtilde*pprime*cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq))
                                          )
                            - (pprime*ktilde)/manderstan_t*sin_theta_kq*sin_theta_pprimeq*f0_E1*f0_E2*(1. + f0_E3)
                             *( (1. - f0_E2)*deltaf_chi(pprime)*(3.*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq*cos_theta_pprimeq)
                               +f0_E3*deltaf_chi(pprime+omega)*(3.*Power(1./(pprime+omega), 2)*(pprime*qtilde*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq + pprime*pprime*sin_theta_pprimeq*cos_theta_pprimeq*sin_theta_kq*cos_theta_kq))
                              );
      equilibrium_result += eq_integrand*pprime_weight[i];
      viscous_result += vis_integrand*pprime_weight[i];
   }
   term_st = equilibrium_result;
   term_st_vis = viscous_result;
   //u/t term
   equilibrium_result = 0.0;
   viscous_result = 0.0;
   for(int i = 0; i < n_pprime; i++)
   {
      double pprime = pprime_pt[i];
      double cos_theta_pprimeq = (omega*omega - qtilde*qtilde + 2.*omega*pprime)/(2.*pprime*qtilde);
      double sin_theta_pprimeq = sqrt(1. - cos_theta_pprimeq*cos_theta_pprimeq);
      double f0_E1 = Fermi_distribution(omega+ktilde);
      double f0_E2 = Fermi_distribution(pprime);
      double f0_E3 = Bose_distribution(omega+pprime);
      double eq_integrand = (- 2.*pprime*ktilde)/manderstan_t*(1. - cos_theta_pprimeq*cos_theta_kq)*f0_E1*f0_E2*(1. + f0_E3);
      double vis_integrand = eq_integrand*((1. - f0_E1)*deltaf_chi(omega+ktilde)*(-0.5 + 1.5*Power((qtilde*cos_theta_kq + ktilde)/(omega + ktilde), 2))
                                          +(1. - f0_E2)*deltaf_chi(pprime)*(-0.5 + 1.5*(cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 0.5*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq))
                                          +f0_E3*deltaf_chi(pprime+omega)*(-0.5 + 1.5*Power(1./(pprime+omega), 2)*(0.5*pprime*pprime*sin_theta_kq*sin_theta_kq*sin_theta_pprimeq*sin_theta_pprimeq + qtilde*qtilde*cos_theta_kq*cos_theta_kq + pprime*pprime*cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq + 2.*qtilde*pprime*cos_theta_kq*cos_theta_kq*cos_theta_pprimeq*cos_theta_pprimeq))
                                          )
                            + (pprime*ktilde)/manderstan_t*sin_theta_kq*sin_theta_pprimeq*f0_E1*f0_E2*(1. + f0_E3)
                             *( (1. - f0_E2)*deltaf_chi(pprime)*(3.*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq*cos_theta_pprimeq)
                               +f0_E3*deltaf_chi(pprime+omega)*(3.*Power(1./(pprime+omega), 2)*(pprime*qtilde*sin_theta_kq*cos_theta_kq*sin_theta_pprimeq + pprime*pprime*sin_theta_pprimeq*cos_theta_pprimeq*sin_theta_kq*cos_theta_kq))
                              );
      equilibrium_result += eq_integrand*pprime_weight[i];
      viscous_result += vis_integrand*pprime_weight[i];
   }
   term_ut = equilibrium_result;
   term_ut_vis = viscous_result;

   results[0] = term_st + term_ut;
   results[1] = term_st_vis + term_ut_vis;
}

void QGP_2to2_Scattering::calculateEmissionrates_soft(string filename_in)
{
   double hbarC = Phycons.get_hbarC();
   filename = filename_in;

   double* results = new double [2];
   double p_min = 0.0;
   double p_max = qtilde_cutoff;
   gauss_quadrature(n_pSoft, 1, 0.0, 0.0, p_min, p_max, pSoft, pSoft_weight);

   for(int i=0; i<n_ktilde; i++)
   {
       double ktilde = ktilde_pt[i];
       double f_q = Fermi_distribution(ktilde);
       double e_sq = Phycons.get_e_sq();
       double q_sq = Phycons.get_q_sq();
       double N_c = Phycons.get_N_c();
       double prefactor = 1./(2.*(8.*M_PI*M_PI*M_PI))
                          *(- e_sq)*q_sq*N_c*8.*f_q/ktilde;
       double equilibrium_result_p = 0.0;
       double viscous_result_p = 0.0;
       for(int k=0; k<n_pSoft; k++)
       {
          double equilibrium_result_theta = 0.0;
          double viscous_result_theta = 0.0;
          for(int l=0; l<n_theta; l++)
          {
             getIntegrand(ktilde, pSoft[k], costhetaSoft[l], results);
             equilibrium_result_theta += results[0]*costhetaSoft_weight[l];
             viscous_result_theta += results[1]*costhetaSoft_weight[l];
          }
          equilibrium_result_p += equilibrium_result_theta*pSoft_weight[k];
          viscous_result_p += viscous_result_theta*pSoft_weight[k];
       }
       equilibriumTilde_results[i] = equilibrium_result_p*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
       viscousTilde_results[i] = ((1. - f_q)*deltaf_chi(ktilde)*equilibrium_result_p + viscous_result_p)*prefactor/(ktilde*ktilde)/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates

   }
   
   buildupEmissionrate2DTable();
   output_emissionrateTable();
   delete [] results;
}

void QGP_2to2_Scattering::getIntegrand(double ktilde, double ptilde, double costhetap, double* results)
{
    double common_factor = ptilde*ptilde/(4.*M_PI*M_PI);
    double equilibrium_integrand;
    double vis_integrand;

    Selfenergy_coefficients Sigma_p;
    Selfenergy_coefficients* Sigma_p_ptr = &Sigma_p;
    double p0tilde = ptilde*costhetap;
    get_quark_selfenergy_coefficients(p0tilde, ptilde, Sigma_p_ptr);
    double ReA0 = Sigma_p.Re_A0;
    double ImA0 = Sigma_p.Im_A0;
    double ReB0 = Sigma_p.Re_B0;
    double ImB0 = Sigma_p.Im_B0;
    double ReA1 = Sigma_p.Re_A1;
    double ImA1 = Sigma_p.Im_A1;
    double ReB1 = Sigma_p.Re_B1;
    double ImB1 = Sigma_p.Im_B1;
    double ReC1 = Sigma_p.Re_C1;
    double ImC1 = Sigma_p.Im_C1;
    //ReA1 = 0.0;
    //ImA1 = 0.0;
    //ReB1 = 0.0;
    //ImB1 = 0.0;
    //ReC1 = 0.0;
    //ImC1 = 0.0;

    //equilibrium part
    double Re_Q_dot_peq0 = -(ktilde*ReB0);
    double Im_Q_dot_peq0 = -(ktilde*ImB0);
    double Re_peq0_dot_peq0 = - Power(ImB0,2) - 2*costhetap*ImA0*ImB0*ptilde
                      - (-1 + Power(costhetap,2))*Power(ptilde,2)*(Power(ImA0,2) 
                      - Power(-1 + ReA0,2)) + 2*costhetap*ptilde*(-1 + ReA0)*ReB0 
                      + Power(ReB0,2);
    double Im_peq0_dot_peq0 = 2*(-(ImA0*Power(ptilde,2)*(-1 + ReA0)) 
                      + Power(costhetap,2)*ImA0*Power(ptilde,2)*(-1 + ReA0) 
                      + ImB0*ReB0 + costhetap*ptilde*(ImB0*(-1 + ReA0) 
                      + ImA0*ReB0));
    double Im_Eqpart = Impart_ComplexDivide(Re_Q_dot_peq0, Im_Q_dot_peq0, 
                                            Re_peq0_dot_peq0, Im_peq0_dot_peq0);
    double Re_Eqpart = Repart_ComplexDivide(Re_Q_dot_peq0, Im_Q_dot_peq0, 
                                            Re_peq0_dot_peq0, Im_peq0_dot_peq0);
    equilibrium_integrand = Im_Eqpart;

    //viscous part
    double Re_Q_dot_Sigma1 = (ktilde*ptilde*((-1 + 3*Power(costhetap,2))*ptilde*ReB1 
                             + 4*costhetap*ReC1))/2.;
    double Im_Q_dot_Sigma1 = (ktilde*ptilde*(4*costhetap*ImC1 - ImB1*ptilde
                             + 3*Power(costhetap,2)*ImB1*ptilde))/2.;
    double Re_peq0_dot_Sigma1 = ((-1 + 3*Power(costhetap,2))*Power(ptilde,2)*(ImB0*(ImB1 
                 + costhetap*ImA1*ptilde) + ImA0*(2*ImC1 + ptilde*(costhetap*ImB1 
                 - ImA1*ptilde + Power(costhetap,2)*ImA1*ptilde)) - Power(ptilde,2)*ReA1 
                 + Power(costhetap,2)*Power(ptilde,2)*ReA1 + Power(ptilde,2)*ReA0*ReA1 
                 - Power(costhetap,2)*Power(ptilde,2)*ReA0*ReA1 
                 - costhetap*ptilde*ReA1*ReB0 + costhetap*ptilde*ReB1 
                 - costhetap*ptilde*ReA0*ReB1 - ReB0*ReB1 + 2*ReC1 
                 - 2*ReA0*ReC1))/2.;
    double Im_peq0_dot_Sigma1 = -((-1 + 3*Power(costhetap,2))*Power(ptilde,2)*(ImA1*Power(ptilde,2) 
                 + 2*ImC1*(-1 + ReA0) - ImA1*Power(ptilde,2)*ReA0 
                 - ImA0*Power(ptilde,2)*ReA1 + Power(costhetap,2)*Power(ptilde,2)
                 *(ImA1*(-1 + ReA0) + ImA0*ReA1) + ImB1*ReB0 + ImB0*ReB1 
                 + costhetap*ptilde*(ImB1*(-1 + ReA0) + ImB0*ReA1 + ImA1*ReB0 
                 + ImA0*ReB1) + 2*ImA0*ReC1))/2.;
    double vis_part1 = - Impart_ComplexDivide(Re_Q_dot_Sigma1, Im_Q_dot_Sigma1, Re_peq0_dot_peq0, Im_peq0_dot_peq0);
    double Re_vis_temp_part2 = Repart_ComplexDivide(Re_peq0_dot_Sigma1, Im_peq0_dot_Sigma1, Re_peq0_dot_peq0, Im_peq0_dot_peq0);
    double Im_vis_temp_part2 = Impart_ComplexDivide(Re_peq0_dot_Sigma1, Im_peq0_dot_Sigma1, Re_peq0_dot_peq0, Im_peq0_dot_peq0);
    double vis_part2 = 2.*Impart_ComplexMultiply(Re_Eqpart, Im_Eqpart, Re_vis_temp_part2, Im_vis_temp_part2);
    vis_integrand = vis_part1 + vis_part2;

    results[0] = common_factor*equilibrium_integrand;
    results[1] = common_factor*vis_integrand;
    return;
}

void QGP_2to2_Scattering::get_quark_selfenergy_coefficients(double p_0_tilde, double p_i_tilde, Selfenergy_coefficients* Sigma_ptr)
{
    double eps = 1e-10;
    double g_s = Phycons.get_g_s_const();
    double C_F = Phycons.get_C_F();

    double omega_0_sq = g_s*g_s*C_F/8.0;

    double p0_tilde_sq = p_0_tilde*p_0_tilde;
    double p_tilde_sq = p_i_tilde*p_i_tilde;
    double p0_tilde_cubic = p0_tilde_sq*p_0_tilde;
    double p_tilde_cubic = p_tilde_sq*p_i_tilde;
    double Re_Qx = quark_selfenergy_Q(p_0_tilde/p_i_tilde);
    double Im_Qx = - 0.5*(p_0_tilde/p_i_tilde)*M_PI;
    
    //Equilibrium quark self energy from Hard Thermal Loop Approximation
    Sigma_ptr->Re_A0 = omega_0_sq/p_tilde_sq*(Re_Qx - 1.); 
    Sigma_ptr->Re_B0 = omega_0_sq*((- p_0_tilde/p_tilde_sq + 1./p_0_tilde)*Re_Qx + p_0_tilde/p_tilde_sq);

    double neq_coeff;
    if(fabs(deltaf_alpha - 2) < eps)
       neq_coeff = (4./(M_PI*M_PI))*12.62159748;  // for deltaf_alpha = 2
    else if(fabs(deltaf_alpha - 1) < eps)
       neq_coeff = (4./(M_PI*M_PI))*4.93480220;   // for deltaf_alpha = 1

    //Off-equilibrium corrections coefficients from Hard Loop Approximation
    Sigma_ptr->Re_A1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                  *((5.*p0_tilde_sq - 3.*p_tilde_sq)*Re_Qx - 5.*p0_tilde_sq + 4./3.*p_tilde_sq);
    Sigma_ptr->Re_B1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                  *((- 5.*p0_tilde_cubic + 6.*p_0_tilde*p_tilde_sq - p_tilde_cubic*p_i_tilde/p_0_tilde)*Re_Qx
                     + 5.*p0_tilde_cubic - 13./3.*p_0_tilde*p_tilde_sq);
    Sigma_ptr->Re_C1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_sq*p_tilde_sq)
                  *((p0_tilde_sq - p_tilde_sq)*Re_Qx - p0_tilde_sq + 2./3.*p_tilde_sq);
    
    //imaginary part
    if(p_0_tilde > p_i_tilde || p_0_tilde < - p_i_tilde)
    {
       Sigma_ptr->Im_A0 = 0.0e0;
       Sigma_ptr->Im_B0 = 0.0e0;
       Sigma_ptr->Im_A1 = 0.0e0;
       Sigma_ptr->Im_B1 = 0.0e0;
       Sigma_ptr->Im_C1 = 0.0e0;
    }
    else
    {
       Sigma_ptr->Im_A0 = omega_0_sq/p_tilde_sq*Im_Qx; 
       Sigma_ptr->Im_B0 = omega_0_sq*(- p_0_tilde/p_tilde_sq + 1./p_0_tilde)*Im_Qx;
       Sigma_ptr->Im_A1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                          *(5.*p0_tilde_sq - 3.*p_tilde_sq)*Im_Qx;
       Sigma_ptr->Im_B1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                          *(- 5.*p0_tilde_cubic + 6.*p_0_tilde*p_tilde_sq - p_tilde_cubic*p_i_tilde/p_0_tilde)*Im_Qx;
       Sigma_ptr->Im_C1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_sq*p_tilde_sq)
                          *(p0_tilde_sq - p_tilde_sq)*Im_Qx;
    }
    return;
}

double QGP_2to2_Scattering::quark_selfenergy_Q(double x)
{
    double result; 
    double eps = 1e-100;
    if(fabs(x) > 1.)
      result = 0.5*x*log((x + 1. - eps)/(x - 1. + eps));
    else
      result = 0.5*x*log((x + 1. + eps)/(1. - x + eps));

    return(result);
}

double QGP_2to2_Scattering::Bose_distribution(double Etilde)
{
   return(1.0/(exp(Etilde)-1.0));
}

double QGP_2to2_Scattering::Fermi_distribution(double Etilde)
{
    return(1.0/(exp(Etilde) + 1.0));
}

double QGP_2to2_Scattering::deltaf_chi(double ptilde)
{ 
    return(pow(ptilde, deltaf_alpha));
}

double QGP_2to2_Scattering::Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Re2 - Im1*Im2);
}

double QGP_2to2_Scattering::Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Im2 + Im1*Re2);
}

double QGP_2to2_Scattering::Repart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
{
    return((Re1*Re2 + Im1*Im2)/(Re2*Re2 + Im2*Im2));
}

double QGP_2to2_Scattering::Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
{
    return((-Re1*Im2 + Im1*Re2)/(Re2*Re2 + Im2*Im2));
}

inline double Power(double x, int a)
{
    double result = 1.0;
    for(int i=0; i<a; i++)
       result *= x;
    return(result);
}
