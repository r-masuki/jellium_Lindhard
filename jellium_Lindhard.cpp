// Calculate Lindhard function of HEG(homogeneous electron gas) or 1,2,3 dimensions.

#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <complex>

// complex unit
const std::complex<double> i_cpx = std::complex<double>(0.0, 1.0);

// step function
double step_fcn(double x){
  if(x > 0.0){
    return 1.0;
  }
  else if(x == 0.0){
    return 0.5;
  }
  else{
    return 0.0;
  }
}

// Lindhard function
std::complex<double> chi_Lindhard_1d(double q, double w){
  std::complex<double> tmp1, tmp2, tmp3;
  double q2 = q*q;

  tmp1 = std::log(std::abs( ((w+2.0*q+q2)*(w-2.0*q-q2)) / ((w+2.0*q-q2)*(w-2.0*q+q2)) ));
  tmp2 = i_cpx*M_PI*step_fcn(-w+2.0*q+q2)*step_fcn(w-2.0*q+q2);
  tmp3 = -i_cpx*M_PI*step_fcn(w+2.0*q+q2)*step_fcn(-w-2.0*q+q2);

  return 1.0/(2.0*q) * (tmp1 + tmp2 + tmp3);
}

int main(){
  
  // log files
  std::ofstream flog, fres;

  flog.open("result/jellium_Lindhard.log");

  //-----------------------------------------------------
  // 1. set parameter
  //-----------------------------------------------------

  std::cerr << "1. set parameters." << std::endl;
  flog << "1. set parameters." << std::endl;

  // calculation parameters
  int dimension = 1;
  double q_max = 5.0;
  int N_q = 100;
  double w_max = 3.0;
  int N_w = 100;

  flog << "dimension = " << dimension << std::endl;
  flog << "q_max = " << q_max << std::endl;
  flog << "N_q = " << N_q << std::endl;
  flog << "w_max = " << w_max << std::endl;
  flog << "N_w = " << N_w << std::endl;
  flog << std::endl;

  //-----------------------------------------------------
  // 2. prepare q-mesh and w-mesh
  //-----------------------------------------------------

  std::cerr << "2. prepare q-mesh and w-mesh." << std::endl;
  flog << "2. prepare q-mesh and w-mesh." << std::endl;

  std::vector<double> q(N_q+1), w(N_w);
  for(int i = 0; i < N_q; i++){
    q[i] = q_max/N_q * (i+1);
  }
  for(int i = 0; i < N_w; i++){
    w[i] = w_max/N_w * i;
  }

  flog << std::endl;

  //-----------------------------------------------------
  // 3. calculate Lindhard function
  //-----------------------------------------------------

  std::cerr << "3. calculate Lindhard function." << std::endl;
  flog << "3. calculate Lindhard function." << std::endl;

  std::vector<std::vector<std::complex<double> > > chi_Lindhard(N_q, std::vector<std::complex<double> >(N_w));

  // 1d 
  if(dimension == 1){
    for(int i = 0; i < N_q; i++){
      for(int j = 0; j < N_w; j++){
        chi_Lindhard[i][j] = chi_Lindhard_1d(q[i], w[j]);
      }
    }
  }

  flog << std::endl;

  //-----------------------------------------------------
  // 4. print result.
  //-----------------------------------------------------
  std::cerr << "4. print result." << std::endl;
  flog << "4. print result." << std::endl;
  
  // print result of chi(q,w=0)
  fres.open("result/q_chi_w=0.txt");
  fres << "# q chi_Lindhard" << std::endl;
  for(int i = 0; i < N_q; i++){
    fres << q[i] << " " << chi_Lindhard[i][0].real() << " " << chi_Lindhard[i][0].imag() << std::endl;
  }
  fres.close();

  // print result of chi(q, w)
  fres.open("result/q_w_chi.txt");
  fres << "# q w chi_Lindhard" << std::endl;
  for(int i = 0; i < N_q; i++){
    for(int j = 0; j < N_w; j++){
      fres << q[i] << " " << w[j] << " " << chi_Lindhard[i][j].real() << " " << chi_Lindhard[i][j].imag() << std::endl;
    }
    fres << std::endl;
  }
  fres.close();
  
  std::cerr << "done.\nend program." << std::endl;
  flog << "done.\\end program." << std::endl;

  flog.close();

  return 0;

}
