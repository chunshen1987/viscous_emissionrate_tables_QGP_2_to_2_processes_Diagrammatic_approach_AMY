#ifndef QGP_2to2_Scattering_H
#define QGP_2to2_Scattering_H

#include <string>
#include "ParameterReader.h"
#include "Physicalconstants.h"

typedef struct
{
   double Re_A0;
   double Im_A0;
   double Re_B0;
   double Im_B0;
   double Re_A1;
   double Im_A1;
   double Re_B1;
   double Im_B1;
   double Re_C1;
   double Im_C1;
}Selfenergy_coefficients;

inline double Power(double x, int a);

class QGP_2to2_Scattering
{
   private:
      ParameterReader *paraRdr;

      Physicalconstants Phycons;

      int n_Eq, n_Temp;
      double *Eq_tb, *T_tb;
      double **equilibrium_results, **viscous_results;

      int channel;
      string filename;
      double qtilde_cutoff;

      int n_ktilde;
      double *ktilde_pt;
      double *equilibriumTilde_results, *viscousTilde_results;

      // Gaussian quadrature points for phase space integrations 
      int n_qtilde_I1;
      double *qtilde_I1_pt, *qtilde_I1_weight, *qtilde_I1_pt_standard, *qtilde_I1_weight_standard;
      int n_qtilde_I2_1;
      double *qtilde_I2_1_pt, *qtilde_I2_1_weight, *qtilde_I2_1_pt_standard, *qtilde_I2_1_weight_standard;
      int n_qtilde_I2_2;
      double *qtilde_I2_2_pt, *qtilde_I2_2_weight, *qtilde_I2_2_pt_standard, *qtilde_I2_2_weight_standard;

      int n_omega_I1;
      double *omega_I1_pt, *omega_I1_weight, *omega_I1_pt_standard, *omega_I1_weight_standard;
      int n_omega_I2;
      double *omega_I2_pt, *omega_I2_weight, *omega_I2_pt_standard, *omega_I2_weight_standard;

      int n_pprime_I1;
      double pprime_max_I1;
      double *pprime_I1_pt, *pprime_I1_weight, *pprime_I1_pt_standard, *pprime_I1_weight_standard;
      int n_pprime_I2;
      double *pprime_I2_pt, *pprime_I2_weight, *pprime_I2_pt_standard, *pprime_I2_weight_standard;

      // Gaussian quadrature points for soft momentum region
      int n_pSoft;
      double *pSoft, *pSoft_weight;
      int n_theta;
      double *costhetaSoft, *costhetaSoft_weight;

      double deltaf_alpha;

   public:
      QGP_2to2_Scattering(ParameterReader* paraRdr_in);
      ~QGP_2to2_Scattering();

      void buildupEmissionrate2DTable();
      void output_emissionrateTable();

      void set_gausspoints();
      
      void calculateEmissionrates_I1(int channel_in, string filename_in);
      void scale_gausspoints_omega(double ktilde);
      void Integrate_I1_qtilde(double ktilde, double omega, double* results);
      void Integrate_I1_pprime(double ktilde, double omega, double qtilde, double* results);

      void calculateEmissionrates_I2(int channel, string filename);
      void scale_gausspoints_qtilde(double ktilde);
      void Integrate_I2_omega(double ktilde, double qtilde, double* results);
      void Integrate_I2_pprime(double ktilde, double qtilde, double omega, double* results);

      void calculateEmissionrates_soft(string filename);
      void getIntegrand(double qtilde, double ptilde, double costhetap, double* results);
      void get_quark_selfenergy_coefficients(double p_0_tilde, double p_i_tilde, Selfenergy_coefficients* Sigma_ptr);
      double quark_selfenergy_Q(double x);

      double Bose_distribution(double Etilde);
      double Fermi_distribution(double Etilde);
      double deltaf_chi(double ptilde);

      double Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
      double Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
      double Repart_ComplexDivide(double Re1, double Im1, double Re2, double Im2);
      double Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2);
};

#endif
