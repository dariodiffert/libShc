//
// Copyright (C) 2017 Dario Differt
// This file is part of libShc (spherical harmonics computations)
// 
// libShc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// libShc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with libShc.  If not, see <http://www.gnu.org/licenses/>.
//

#include <Shc.h>

using namespace std;
using namespace Eigen;
using namespace shc;
 
float Shc::feature_difference(VectorReal& v1, VectorReal& v2) {
  
  float n0=0;
  VectorReal v0;   
  
  if (v1.size() != v2.size()) {
    print_warning("feature_difference", "This should never be reached, bug?");
  }
    
  switch (feature_norm) {
    case NORM_1:
      v0 = v1-v2;
      n0  = v0.lpNorm<1>();
      break;
    case NORM_2:
      v0 = v1-v2;
      n0  = v0.lpNorm<2>();
      break;
    case NORM_INFINITY:
      v0 = v1-v2;
      n0  = v0.lpNorm<Eigen::Infinity>();
      break;
  }
  
  if (bool_feature_normalize == false) {
    
    return n0;
    
  } else {
    
    float n1=0; float n2=0;  
    switch (feature_norm) {
      case NORM_1:
        n1  = v1.lpNorm<1>();
        n2  = v2.lpNorm<1>();
        break;
      case NORM_2:
        n1  = v1.lpNorm<2>();
        n2  = v2.lpNorm<2>();
        break;
      case NORM_INFINITY:
        n1  = v1.lpNorm<Eigen::Infinity>();
        n2  = v2.lpNorm<Eigen::Infinity>();
        break;
    }
    
    if (n1+n2 == 0) {
      return 1; // return max difference if both are 0
    } else {
      return n0 / (n1+n2);
    }
    
  }
  
}
float Shc::calc_difference_ISE(Shpm& shpm1, Shpm& shpm2) {
  
  Coef coef1 = calc_feature_ISE(shpm1);
  Coef coef2 = calc_feature_ISE(shpm2);
  
  return feature_difference(coef1, coef2);
  
}
float Shc::calc_difference_AS(Shpm& shpm1, Shpm& shpm2) {
  
  Coef coef1 = calc_feature_amplitudespectrum(shpm1);
  Coef coef2 = calc_feature_amplitudespectrum(shpm2);
  
  return feature_difference(coef1, coef2);
  
}
float Shc::calc_difference_BS(Shpm& shpm1, Shpm& shpm2) {
  
  Coef coef1 = calc_feature_bispectrum(shpm1);
  Coef coef2 = calc_feature_bispectrum(shpm2);
  
  return feature_difference(coef1, coef2);
  
}
float Shc::calc_difference_BS_dense(Shpm& shpm1, Shpm& shpm2) {
  
  Coef coef1 = calc_feature_bispectrum_dense(shpm1);
  Coef coef2 = calc_feature_bispectrum_dense(shpm2);
  
  return feature_difference(coef1, coef2);
  
}

VectorReal Shc::calc_feature_ISE(Shpm& shpm) {
  
  int i_size  = l2sh(shpm.l_max, shpm.contin);
  
  VectorReal ISE = shpm.coef.segment(0, i_size);
  
  return shpm.coef;
  
}
VectorReal Shc::calc_feature_amplitudespectrum(Shpm& shpm) {
  
  float sum = 0;
  int countRes = 0;
  
  VectorReal v(shpm.l_max);

  for (int l=0; l<shpm.l_max; l++) {
    
    if (contin_skip(shpm.contin, l) == true) {continue;}

    int i = l2sh(l, shpm.contin);
    
    sum = 0;
    for (int j=i; j<i+2*l+1; j++) {
      sum += shpm.coef[j]*shpm.coef[j];
    }
    
    v(countRes) = sum;
    countRes++;
    
  } 
  
  return v.segment(0,countRes);
  
}
VectorReal Shc::calc_feature_bispectrum(Shpm& shpm) {
  
  int l_max = min(n_bands_CG, shpm.l_max); 
  int n = l_max*l_max + (l_max-1)*(2*l_max-1)*l_max / 6;
  
  VectorReal result(n); // reserve maximal desired space (its always required less)
  result.setZero();
  
  CoefBandwise bandwise = shpm2bandwise(shpm, l_max);
  
  // calculate bispectrum
  int countRes = 0;
  for (int l1=0; l1<l_max; l1++) {
    for (int l2=l1; l2<l_max; l2++) {
      
      if (contin_skip(shpm.contin, l1) || contin_skip(shpm.contin, l2)) {continue;}
      
      // calculate [Fk tensor Fl]
      VectorReal Fkl((2*l1+1)*(2*l2+1));
      Fkl = kroneckerProduct(bandwise[l1],bandwise[l2]);

      // calculate bispectrum for each i
      for (int i=l2-l1; i<=l2+l1; i++) {
        
        if (contin_skip(shpm.contin, i)) {continue;}
        
        if (i < l_max) {
          int np   = i*i - (l1-l2)*(l1-l2);
          int cols = 2*i+1;
          int rows = lut_cg_realCoupling[l1][l2].rows();
          result(countRes) = bandwise[i].transpose() * lut_cg_realCoupling[l1][l2].block(0,np,rows,cols).transpose() * Fkl;
        } else {
          result(countRes) = 0;
        }
        countRes++;
        
      }

    }
  }
  
  return result.segment(0, countRes);
   
}
VectorReal Shc::calc_feature_bispectrum_dense(Shpm& shpm) {
  
  if (trifold(shpm.contin) != FULL) {
    print_warning("calc_feature_bispectrum_dense", "Only Shpm without hemispherical continuation (i.e. contin == FULL) can be used.");
    return VectorReal();
  }
  
  int l_max = min(n_bands_CG, shpm.l_max); 
  int n = l_max*l_max + (l_max-1)*(2*l_max-1)*l_max / 6;
  
  VectorReal result(n); // reserve maximal desired space (its always required less)
  result.setZero();
  
  CoefBandwise bandwise = shpm2bandwise(shpm, l_max);
  
  // calculate bispectrum
  int countRes = 0; int l2; int i;
  VectorReal Fkl;
  for (int l1=1; l1<l_max; l1++) {
      
    l2 = l1-1;
    i = 1;
    Fkl((2*l1+1)*(2*l2+1));
    Fkl = kroneckerProduct(bandwise[l1],bandwise[l2]);
    if (i < l_max) {
      int np   = i*i - (l1-l2)*(l1-l2);
      int cols = 2*i+1;
      int rows = lut_cg_realCoupling[l1][l2].rows();        
      result(countRes) = bandwise[i].transpose() * lut_cg_realCoupling[l1][l2].block(0,np,rows,cols).transpose() * Fkl;
    } else {
      result(countRes) = 0;
    }
    countRes++;

  }
  
  return result.segment(0, countRes);
   
}


