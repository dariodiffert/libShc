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

Shc::s_compass_prepare_shpm Shc::rotate_prepare_shpm(VecShpm& vshpm) {
  
  s_compass_prepare_shpm scp;
  scp.n_mask = 0;
  scp.n_shpm = vshpm.size();
  scp.n_total = scp.n_shpm;
  scp.cv.reserve(scp.n_total);
  for (int i=0; i<scp.n_shpm; i++) {
    scp.cv.push_back(vshpm[i]);
  }
  scp.cv_cur.resize(scp.n_total);
  
  scp.l_max.resize(scp.n_total);  
  scp.sh_max.resize(scp.n_total);  
  scp.valid.resize(scp.n_total);
  scp.rr_cur.resize(scp.n_total);
  
  for (int q=0; q<scp.n_total; q++) {
    scp.valid[q] = true;
  }
  
  for (int q=0; q<scp.n_total; q++) {
    scp.l_max[q].resize(n_rot_par);
    scp.sh_max[q].resize(n_rot_par);
    for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
      scp.l_max[q][index_rot_par]  = min(scp.cv[q].l_max, v_rotation[index_rot_par].l_max);
      scp.sh_max[q][index_rot_par] = l2sh(scp.l_max[q][index_rot_par], scp.cv[q].contin);
    }
  }
  
  return scp;
  
}
Shc::s_compass_prepare_shpm Shc::rotate_prepare_shpm(Shpm& shpm) {
  
  s_compass_prepare_shpm scp;
  scp.n_mask = 0;
  scp.n_shpm = 1;
  scp.n_total = 1;
  scp.cv.push_back(shpm);
  scp.cv_cur.resize(scp.n_total);
  
  scp.l_max.resize(scp.n_total);  
  scp.sh_max.resize(scp.n_total);  
  scp.valid.resize(scp.n_total);
  scp.rr_cur.resize(scp.n_total);
  
  scp.valid[0] = true;
  
  scp.l_max[0].resize(n_rot_par);
  scp.sh_max[0].resize(n_rot_par);
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    scp.l_max[0][index_rot_par]  = min(scp.cv[0].l_max, v_rotation[index_rot_par].l_max);
    scp.sh_max[0][index_rot_par] = l2sh(scp.l_max[0][index_rot_par], scp.cv[0].contin);
  }
  
  return scp;
  
}
Shc::s_compass_prepare_shpm Shc::compass_prepare_shpm(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2) {
  
  s_compass_prepare_shpm scp;
  scp.n_mask = v_mask1.size();
  scp.n_shpm = v_shpm1.size();
  scp.n_total = scp.n_mask + scp.n_shpm;
  scp.cv.reserve(scp.n_total);
  scp.ss.reserve(scp.n_total);
  for (int i=0; i<scp.n_shpm; i++) {
    scp.cv.push_back(v_shpm1[i]);
    scp.ss.push_back(v_shpm2[i]);
  }
  for (int i=0; i<scp.n_mask; i++) {
    scp.cv.push_back(v_mask1[i]);
    scp.ss.push_back(v_mask2[i]);
  }
  scp.cv_cur.resize(scp.n_total);
  scp.ss_cur.resize(scp.n_total);
  
  scp.l_max.resize(scp.n_total);  
  scp.sh_max.resize(scp.n_total);  
  scp.valid.resize(scp.n_total);
  scp.rr_cur.resize(scp.n_total);
  
  for (int q=0; q<scp.n_total; q++) {
    scp.valid[q] = true;
    if ((scp.cv[q].contin == HEMI_RM && scp.ss[q].contin == HEMI_RMN) || (scp.cv[q].contin == HEMI_RMN && scp.ss[q].contin == HEMI_RM)) {
      scp.valid[q] = false;
    } else if (trifold(scp.cv[q].contin) != FULL && trifold(scp.ss[q].contin) == FULL) {
      scp.ss[q] = contin_convert(scp.ss[q], scp.cv[q].contin);
    } else if (trifold(scp.ss[q].contin) != FULL && trifold(scp.cv[q].contin) == FULL) {
      scp.cv[q] = contin_convert(scp.cv[q], scp.ss[q].contin);
    } 
  }
  
  int l_max_rotPar_value = 0;
  int l_max_rotPar_index = 0;
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (v_rotation[index_rot_par].l_max_compass > l_max_rotPar_value) {
      l_max_rotPar_value = v_rotation[index_rot_par].l_max_compass;
      l_max_rotPar_index = index_rot_par;
    }
  }
  
  for (int q=0; q<scp.n_total; q++) {
    scp.l_max[q].resize(n_rot_par);
    scp.sh_max[q].resize(n_rot_par);
    for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
      scp.l_max[q][index_rot_par]  = min3(scp.cv[q].l_max, scp.ss[q].l_max, v_rotation[index_rot_par].l_max_compass);      
      scp.sh_max[q][index_rot_par] = l2sh(scp.l_max[q][index_rot_par], scp.cv[q].contin);
    }
  }
    
  // instead of calculating resultZXZXYZY we only calculate resultZXZXYZ0 and rotate the ss once instead
  compass_calc_SS(scp.ss, scp, l_max_rotPar_index);
  
  switch (tangent_distance_mode) {
    case OFF:
      break;
    case ONESIDED:
      scp.ss_td.resize(scp.n_shpm);
      scp.L.resize(scp.n_shpm);
      scp.LL.resize(scp.n_shpm);
      for (int i=0; i<scp.n_shpm; i++) {
        scp.ss_td[i].resize(n_rot_par);
        scp.L[i].resize(n_rot_par);
        scp.LL[i].resize(n_rot_par);
        for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
          scp.ss_td[i][index_rot_par].contin = scp.ss[i].contin;
          scp.ss_td[i][index_rot_par].l_max  = scp.l_max[i][index_rot_par];
          scp.ss_td[i][index_rot_par].coef   = scp.ss[i].coef.segment(0, scp.sh_max[i][index_rot_par]);
          tangent_distance_prepare(scp.ss_td[i][index_rot_par], scp.L[i][index_rot_par], scp.LL[i][index_rot_par], v_rotation[index_rot_par].td_active);
        }
      }
      break;
    case TWOSIDED:
      print_warning("compass_prepare_shpm", "Compass cannot be used with tangent_distance mode TWOSIDED.");
      break;
  }
  

  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    scp.rr_cur[q].contin = scp.cv[q].contin;
    scp.rr_cur[q].l_max  = scp.cv[q].l_max;
    scp.rr_cur[q].coef.resize(scp.cv[q].coef.size());
  }
  
  return scp;
  
}
Xyz Shc::intern_compass(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2) {
  
  vector<Xyz> xyz_results;
  Xyz xyz(0,0,0);
  bool rotations_performed = false;
  
  s_compass_prepare_shpm scp = compass_prepare_shpm(v_shpm1, v_shpm2, v_mask1, v_mask2);
  s_rotation_indices indices;
  
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {

    if (v_rotation[index_rot_par].l_max_compass == 0) {
      xyz = Xyz(0,0,0);
    } else {
      
      for (int q=0; q<scp.n_total; q++) {
        scp.cv_cur[q] = scp.cv[q];
        if (tangent_distance_mode == ONESIDED && q < scp.n_shpm) {
          scp.ss_cur[q] = scp.ss_td[q][index_rot_par];
        } else {
          scp.ss_cur[q] = scp.ss[q];
        }
        
        if (scp.cv_cur[q].l_max > scp.l_max[q][index_rot_par] || scp.cv_cur[q].coef.size() != scp.sh_max[q][index_rot_par]) {
          scp.cv_cur[q].l_max = scp.l_max[q][index_rot_par];
          scp.cv_cur[q].coef  = scp.cv[q].coef.segment(0, scp.sh_max[q][index_rot_par]).eval();
        }
        scp.rr_cur[q] = scp.cv_cur[q];
        if (scp.ss_cur[q].l_max > scp.l_max[q][index_rot_par] || scp.ss_cur[q].coef.size() != scp.sh_max[q][index_rot_par]) {
          scp.ss_cur[q].l_max = scp.l_max[q][index_rot_par];
          scp.ss_cur[q].coef  = scp.ss[q].coef.segment(0, scp.sh_max[q][index_rot_par]).eval();
        }
      }
      
      indices = compass_evaluate(index_rot_par, scp);
      rotations_performed = true;
    }

    xyz = rotationIndices2xyz(index_rot_par, indices);
    xyz_results.push_back(xyz);
    
  }  
  
  MatrixRotation R = MatrixReal::Identity(3,3);
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    MatrixRotation RR = xyz2r(xyz_results[index_rot_par]);
    R = RR * R;
  }
  xyz = r2xyz(R);
  
  if (rotations_performed == false) {
    print_warning("compass", "No rotations for compass use initialized!");
  }
  
  return xyz;
  
} 
void Shc::compass_calc_SS(VecShpm& result, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp = result[q].coef;
    
    int l_max = min(scp.l_max[q][index_rot_par], result[q].l_max);
    
    count = 0;
    for (int l=0; l<l_max; l++) {
      if (contin_skip(result[q].contin, l)) {continue;};
      result[q].coef.segment(count,   l)   = v_rotation[index_rot_par].rotate_YN.M1[l] * temp.segment(count,   l);
      result[q].coef.segment(count+l, l+1) = v_rotation[index_rot_par].rotate_YN.M2[l] * temp.segment(count+l, l+1);
      count += 2*l+1;
    }
    
  }
  
}
void Shc::compass_calc_Z(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);

    count = 0;
    if (v_rotation[index_rot_par].psi[index] == 0) {
      temp = in[q].coef.segment(0, scp.sh_max[q][index_rot_par]);
    } else {
      for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
        if (contin_skip(in[q].contin, l)) {continue;};
        temp.segment(count, 2*l+1) = in[q].coef.segment(count, 2*l+1).array() * v_rotation[index_rot_par].rotate_Z[index].v1[l].array() + in[q].coef.segment(count, 2*l+1).reverse().array() * v_rotation[index_rot_par].rotate_Z[index].v2[l].array();
        count += 2*l+1;
      }
    }

    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      result[q].coef.segment(count, 2*l+1) = v_rotation[index_rot_par].rotate_XP.M[l] * temp.segment(count, 2*l+1);
      count += 2*l+1;
    }
    
  }
  
}
void Shc::compass_calc_Y(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);
    
    count = 0;
    if (v_rotation[index_rot_par].psi[index] == 0) {
      temp = in[q].coef.segment(0, scp.sh_max[q][index_rot_par]);
    } else {
      for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
        if (contin_skip(in[q].contin, l)) {continue;};
        temp.segment(count, 2*l+1) = in[q].coef.segment(count, 2*l+1).array() * v_rotation[index_rot_par].rotate_Z[index].v1[l].array() + in[q].coef.segment(count, 2*l+1).reverse().array() * v_rotation[index_rot_par].rotate_Z[index].v2[l].array();
        count += 2*l+1;
      }
    }
    
    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      result[q].coef.segment(count, 2*l+1) = v_rotation[index_rot_par].rotate_YNXN.M[l] * temp.segment(count, 2*l+1);
      count += 2*l+1;
    }
    
  }
  
}
void Shc::compass_calc_X(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    count = 0;
    if (v_rotation[index_rot_par].psi[index] == 0) {
      result[q].coef.segment(0, scp.sh_max[q][index_rot_par]) = in[q].coef.segment(0, scp.sh_max[q][index_rot_par]);
    } else {
      for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
        if (contin_skip(in[q].contin, l)) {continue;};
        result[q].coef.segment(count, 2*l+1) = in[q].coef.segment(count, 2*l+1).array() * v_rotation[index_rot_par].rotate_Z[index].v1[l].array() + in[q].coef.segment(count, 2*l+1).reverse().array() * v_rotation[index_rot_par].rotate_Z[index].v2[l].array();
        count += 2*l+1;
      }
    }
    
  }
  
}
void Shc::compass_calc_CV(VecShpm& result, VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    int l_max = min(scp.l_max[q][index_rot_par], result[q].l_max);
    
    count = 0;
    for (int l=0; l<l_max; l++) {
      if (contin_skip(result[q].contin, l)) {continue;};
      result[q].coef.segment(count,   l)   = v_rotation[index_rot_par].rotate_YP.M1[l] * in[q].coef.segment(count,   l);
      result[q].coef.segment(count+l, l+1) = v_rotation[index_rot_par].rotate_YP.M2[l] * in[q].coef.segment(count+l, l+1);
      count += 2*l+1;
    }
    
  }
  
}
float Shc::compass_evaluate_compare(s_compass_prepare_shpm& scp, int index_rot_par) {
  
  float sum = 0;
  float cur;
  
  for (int q=0; q<scp.n_shpm; q++) {
    if (!scp.valid[q]) {continue;}
    
    if (tangent_distance_mode == ONESIDED) {
      tangent_distance_calc_coef(scp.ss_cur[q].coef, scp.rr_cur[q], scp.ss_cur[q], scp.L[q][index_rot_par], scp.LL[q][index_rot_par]);
    }
    
    if (scp.n_mask == 0) {
      
      cur = feature_difference(scp.rr_cur[q].coef, scp.ss_cur[q].coef);
      
    } else {
            
      Shpm c  = linear_combination(1, scp.rr_cur[q], -1, scp.ss_cur[q]);  
      Shpm m  = intern_product(scp.rr_cur[q+scp.n_shpm], scp.ss_cur[q+scp.n_shpm]); 
      Shpm cm = intern_product(m, c);
      
      float p;
      float sum = 0;
      int n1 = min(c.coef.size(), cm.coef.size());
      for (int i=0; i<n1; i++) {
        p = c.coef(i) * cm.coef(i); // numerator: weighted difference
        if (p>=0) {sum += p;} // only use positive results of the scalar product, otherwise it can become negative
      }
      
      float sumM = 0;
      int n2 = min(scp.rr_cur[q+scp.n_shpm].coef.size(), scp.ss_cur[q+scp.n_shpm].coef.size());
      for (int i=0; i<n2; i++) {
        p = scp.rr_cur[q+scp.n_shpm].coef(i) * scp.ss_cur[q+scp.n_shpm].coef(i); // denominator: sum of weights (normalization)
        if (p>=0) {sumM += p;}
      }
      
      if (sumM > 0) {
        cur = sum / sumM;
      } else {
        cur = sum;
      }
      
    }
    sum += cur;
  }
  
  return sum;
        
}
Shc::s_rotation_indices Shc::compass_evaluate(int index_rot_par, s_compass_prepare_shpm& scp) {
  
  VecShpm Z = scp.cv_cur;
  VecShpm Y = scp.cv_cur;
  VecShpm bestMatch;
  
  float min_diff = 1e12;
  s_rotation_indices indices;
  
  // calculate compass
  int count = 0;
  for (int i=0; i<v_rotation[index_rot_par].n_fast_indices_z; i++) {
    compass_calc_Z(Z, scp.cv_cur, v_rotation[index_rot_par].fast_indices_z[i], scp, index_rot_par);
    
    for (int j=0; j<v_rotation[index_rot_par].n_fast_indices_y; j++) {
      compass_calc_Y(Y, Z, v_rotation[index_rot_par].fast_indices_y[j], scp, index_rot_par);
      
      for (int k=0; k<v_rotation[index_rot_par].n_fast_indices_x; k++) {
        compass_calc_X(scp.rr_cur, Y, v_rotation[index_rot_par].fast_indices_x[k], scp, index_rot_par);
        
        float diff = compass_evaluate_compare(scp, index_rot_par);
        if (diff < min_diff) {
          min_diff = diff;
          bestMatch = scp.rr_cur;
          indices.i = i;
          indices.j = j;
          indices.k = k;
        }
        
        count++;

      }
    }
  }
  
  compass_calc_CV(scp.cv, bestMatch, scp, index_rot_par);
  
  return indices;
  
}

