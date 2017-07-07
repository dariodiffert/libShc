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

void Shc::tangent_distance_prepare(Shpm& shpm, MatrixReal& L, VectorInt active) {
  
  vector<Shpm> T;

  Shpm temp;
  for (uint j=0; j<(uint)tangent_distance.unit.size(); j++) {
    
    bool proceed = false;
    if (active.size() == 0) {
      proceed = true;
    } else {
      for (uint i=0; i<(uint)active.size(); i++) {
        if ((uint)active(i) == j) {
          proceed = true;
          break;
        }
      }
    }
    
    if (proceed == true) {
      temp = intern_transform(shpm, tangent_distance.unit[j]);
      T.push_back(temp);
    }
  }

  int m = T.size();
  int n = shpm.coef.size();

  L.resize(n,m);
  for (int j=0; j<m; j++) {
    L.block(0,j,n,1) = T[j].coef-shpm.coef;
  }
  
}
void Shc::tangent_distance_prepare(Shpm& shpm, MatrixReal& L, MatrixReal& LL, VectorInt active) {
  
  tangent_distance_prepare(shpm, L, active);
  
  MatrixReal LTL;
  MatrixReal LTLI;

  LTL   = L.transpose() * L;
  LTLI  = LTL.inverse();
  LL = LTLI * L.transpose();
  
}

void Shc::tangent_distance_calc_coef(Coef& result, Shpm& shpmE, Shpm& shpmP, MatrixReal& L, MatrixReal& LL) {
  
//   int l_max = min(shpmE.l_max, shpmP.l_max);
//   int n = l_max*l_max;
//   
//   VectorReal alpha = LL.block(0,0,n,n) * (shpmE.coef.segment(0,n) - shpmP.coef.segment(0,n));
//   result = shpmP.coef.segment(0,n) + L.block(0,0,n,n) * alpha;

  VectorReal alpha = LL * (shpmE.coef - shpmP.coef);
  result = shpmP.coef + L * alpha;
  
}
Shpm Shc::tangent_distance_calc(Shpm& shpmE, Shpm& shpmP) {
  
  MatrixReal L, LL;
  tangent_distance_prepare(shpmP, L, LL, VectorInt());
  
  Shpm result = shpmP;
  tangent_distance_calc_coef(result.coef, shpmE, shpmP, L, LL);
  
  return result;
  
}
VecShpm Shc::tangent_distance_calc(VecShpm& shpmE, Shpm& shpmP) {
  
  MatrixReal L, LL;
  tangent_distance_prepare(shpmP, L, LL, VectorInt());
  
  int n = shpmE.size();
  VecShpm result; result.reserve(n);
  
  for (int i=0; i<n; i++) {
    result.push_back(shpmP);
    tangent_distance_calc_coef(result[i].coef, shpmE[i], shpmP, L, LL);
  }
  
  return result;
  
}
void Shc::tangent_distance_calc(Shpm& shpmE, Shpm& shpmP, Shpm& resultE, Shpm& resultP) {
  
  MatrixReal LE; tangent_distance_prepare(shpmE, LE, VectorInt());
  MatrixReal LP; tangent_distance_prepare(shpmP, LP, VectorInt());
  
  MatrixReal LPE = LP.transpose() * LE;
  MatrixReal LEP = LPE.transpose();
  MatrixReal LPP = LP.transpose() * LP;
  MatrixReal LEE = LE.transpose() * LE;
  MatrixReal LEEI = LEE.inverse();
  MatrixReal LPPI = LPP.inverse();
   
  float k1 = 1.0+tangent_distance_spring_constant;
  float k2 = k1*k1;
  MatrixReal LHP = (LPE * LEEI * LE.transpose() - k1 * LP.transpose()) * (shpmE.coef-shpmP.coef);
  MatrixReal RHP = (LPE * LEEI * LEP - k2 * LPP);
  
  MatrixReal LHE = (LEP * LPPI * LP.transpose() - k1 * LE.transpose()) * (shpmE.coef-shpmP.coef);
  MatrixReal RHE = (k2 * LEE - LEP * LPPI * LPE);
  
  MatrixReal alphaP = RHP.inverse() * LHP;
  MatrixReal alphaE = RHE.inverse() * LHE;
  
  resultP = shpmP;
  resultE = shpmE;
  resultP.coef = shpmP.coef + LP * alphaP;
  resultE.coef = shpmE.coef + LE * alphaE;
   
}

