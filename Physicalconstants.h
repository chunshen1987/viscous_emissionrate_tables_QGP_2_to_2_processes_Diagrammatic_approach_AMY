#ifndef Physicalconstants_H
#define Physicalconstants_H

#include "ParameterReader.h"

class Physicalconstants
{
   private:
      ParameterReader *paraRdr;

      double hbarC;
      double alpha_EM;
      double e_sq;
      double q_sq;
      double N_c;
      double C_F;
      double d_F;
      double N_F;

      double g_s;

   public:
      Physicalconstants(ParameterReader* paraRdr_in);
      ~Physicalconstants();

      double get_hbarC() {return(hbarC);};
      double get_alphaEM() {return(alpha_EM);};
      double get_e_sq() {return(e_sq);};
      double get_q_sq();
      double get_C_F() {return(C_F);};
      double get_d_F() {return(d_F);};
      double get_N_c() {return(N_c);};
      double get_N_F() {return(N_F);};
      double get_g_s_const() {return(g_s);};
      double get_g_s_running(double T);
};

#endif
