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

bool Shc::calc_rotate_isEntryZero_X(int mi, int mj, int l_ninety) {
  
  if (mi < 0  && mj < 0  && odd(mi+mj))  {return true;}
  if (mi >= 0 && mj < 0  && even(mi+mj)) {return true;}
  if (mi < 0  && mj >= 0 && even(mi+mj)) {return true;}
  if (mi >= 0 && mj >= 0 && odd(mi+mj))  {return true;}

  if (l_ninety > -1) {
    if (even(l_ninety+mi) && mj  < 0) {return true;}
    if (odd(l_ninety+mi)  && mj >= 0) {return true;}
  }
  
  return false;
  
}
bool Shc::calc_rotate_isEntryZero_Y(int mi, int mj, int l_ninety) {
  
  if (mi >= 0 && mj < 0) {return true;}
  if (mi < 0 && mj >= 0) {return true;}
  
  if (l_ninety > -1) {
    if (even(l_ninety+mi+mj) && mj  < 0)  {return true;}
    if (odd(l_ninety+mi+mj)  && mj >= 0)  {return true;}
  }
  
  return false;
  
}
bool Shc::calc_rotate_isEntryZero_Z(int mi, int mj, int l_ninety) {
  
  if (abs(mi) != abs(mj)) {return true;}
  
  if (l_ninety > -1) {
    if (mi*mj >= 0 && odd(mi))  {return true;}
    if (mi*mj  < 0 && even(mi)) {return true;}
  }
  
  return false;
  
}
bool Shc::calc_rotate_isEntryZero(e_axis axis, int li, int lj, int mi, int mj, bool ninety) {
  
  if (mi < -li || mi > li) {return true;}
  if (mj < -lj || mj > lj) {return true;}
  if (li != lj)            {return true;}
  
  if (ninety) {
    switch (axis) {
      case AXIS_X:
        return calc_rotate_isEntryZero_X(mi, mj, li);
      case AXIS_Y:
        return calc_rotate_isEntryZero_Y(mi, mj, li);
      case AXIS_Z:
        return calc_rotate_isEntryZero_Z(mi, mj, li);
    }
  } else {
    switch (axis) {
      case AXIS_X:
        return calc_rotate_isEntryZero_X(mi, mj, -1);
      case AXIS_Y:
        return calc_rotate_isEntryZero_Y(mi, mj, -1);
      case AXIS_Z:
        return calc_rotate_isEntryZero_Z(mi, mj, -1);
    }
  }
  
  return false;
  
}
bool Shc::calc_rotate_isEntryZero(int li, int lj, int mi, int mj) {
  
  if (mi < -li || mi > li) {return true;}
  if (mj < -lj || mj > lj) {return true;}
  if (li != lj)            {return true;}
    
  return false;
  
}
bool Shc::calc_translate_isEntryZero(int mi, int mj) {
  
  if (mi != mj) {
    return true;
  } else {
    return false;
  }
  
}

MatrixReal Shc::calc_rotate_simple(Xyz xyz, int l_max) {
  
  xyz.x *= -1;
  xyz.y *= -1;
  
  if (l_max < 0) {
    l_max = n_bands;
  }
    
  MatrixRotation RR = xyz2r(xyz);
  MatrixReal M(l2sh(n_bands, FULL), l2sh(n_bands, FULL));
  M.setZero();
  
  MatrixRotation R;

  R(0,0) = RR(1,1); R(0,1) = RR(1,2); R(0,2) = RR(1,0);
  R(1,0) = RR(2,1); R(1,1) = RR(2,2); R(1,2) = RR(2,0);
  R(2,0) = RR(0,1); R(2,1) = RR(0,2); R(2,2) = RR(0,0); 

  int li; int lj; int mi; int mj;
  
  M(0,0) = 1;
  if (n_bands > 1) {
    for (int i=1; i<4; i++) {
      for (int j=1; j<4; j++) {
        index2sh(i, li, mi);
        index2sh(j, lj, mj);
        M(i,j) = R(i-1,j-1);
      }
    }
  }

  int ij_max = l2sh(n_bands, FULL);
  for (int i=4; i<ij_max; i++) {
    for (int j=4; j<ij_max; j++) {
      
      index2sh(i, li, mi);
      index2sh(j, lj, mj);
      
      if (li >= l_max || lj >= l_max) {
        continue;
      }
      
      if (li < 2 || lj < 2) {continue;}
      
      if (calc_rotate_isEntryZero(li, lj, mi, mj)) {continue;}
      
      int x = li*li+li;
      M(x+mi,x+mj) = calc_rotation_entry(R, M, li, mi, mj);
      
    }
  }

  return M;
   
}
MatrixReal Shc::calc_rotate_axis(float rotation_angle, e_axis axis, int l_max) {
  
  if (l_max < 0) {
    l_max = n_bands;
  }
  
  Xyz xyz(0,0,0);
  
  // Internal coordinate system differs from paper
  switch (axis) {
    case AXIS_Z: 
      xyz.z = rotation_angle;
      break;
    case AXIS_Y:
      xyz.y = -rotation_angle;
      break;
    case AXIS_X:
      xyz.x = -rotation_angle;
      break;
    default:
      print("[Shc::rotate] invalid axis label\n");
      break;
  }
  
  bool bool_ninety = false;
  if (angular_difference(rotation_angle, M_PI/2) < 1e-5 || angular_difference(rotation_angle, -M_PI/2) < 1e-5) {
    bool_ninety = true;
  }
  
  MatrixRotation RR = xyz2r(xyz);
  MatrixReal M(l2sh(n_bands, FULL), l2sh(n_bands, FULL));
  M.setZero();
  
  MatrixRotation R;

  R(0,0) = RR(1,1); R(0,1) = RR(1,2); R(0,2) = RR(1,0);
  R(1,0) = RR(2,1); R(1,1) = RR(2,2); R(1,2) = RR(2,0);
  R(2,0) = RR(0,1); R(2,1) = RR(0,2); R(2,2) = RR(0,0); 

  int li; int lj; int mi; int mj;
  
  M(0,0) = 1;
  if (n_bands > 1) {
    for (int i=1; i<4; i++) {
      for (int j=1; j<4; j++) {
        index2sh(i, li, mi);
        index2sh(j, lj, mj);
        if (calc_rotate_isEntryZero(axis, li, lj, mi, mj, bool_ninety)) {continue;}
        M(i,j) = R(i-1,j-1);
      }
    }
  }

  int ij_max = l2sh(n_bands, FULL);
  for (int i=4; i<ij_max; i++) {
    for (int j=4; j<ij_max; j++) {
      
      index2sh(i, li, mi);
      index2sh(j, lj, mj);
      
      if (li >= l_max) {
        break;
      }
      
      if (li < 2 || lj < 2) {continue;}
      if (calc_rotate_isEntryZero(axis, li, lj, mi, mj, bool_ninety)) {continue;}
      
      int x = li*li+li;
      float entry;
      if (mj >= mi) {
        if (axis == AXIS_Y || axis == AXIS_X) {
          entry = calc_rotation_entry(R, M, li, mi, mj);
        } else {
          entry = calc_rotation_entry_Z(xyz.z, mi, mj);
        }
        
        M(x+mi,x+mj) = entry;
        
      } else {
        
        // derived rotation matrix symmetries, nearly halves calculation duration
        if (axis == AXIS_Z) {
          M(x+mi,x+mj) = -M(x+mj,x+mi);
        } else {
          if ((mi+mj)%2 == 0) {
            M(x+mi,x+mj) =  M(x+mj,x+mi);
          } else {
            M(x+mi,x+mj) = -M(x+mj,x+mi);
          }
        }
        
      }

    }
  }

  return M;
   
}
Shc::s_rot_Z Shc::create_rotate_Z(float rotation_angle, int l_max) {
  
  s_rot_Z rot;
  
  rot.v1.resize(l_max);
  rot.v2.resize(l_max);
  
  for (int l=0; l<l_max; l++) {
    
    rot.v1[l].resize(2*l+1);
    rot.v2[l].resize(2*l+1);
    
    for (int m=-l; m<=l; m++) {
      rot.v1[l](m+l) = cos(abs(m)*rotation_angle);
      if (m == 0) {
        rot.v2[l](m+l) = 0;
      } else {
        rot.v2[l](m+l) = sin(-m*rotation_angle);
      }
    }
    
  }

  return rot;
  
}
Shc::s_rot_Y Shc::create_rotate_Y(MatrixReal& M, int l_max) {
  
  s_rot_Y rot;
  
  rot.M1.resize(l_max);
  rot.M2.resize(l_max);
  
  for (int l=0; l<l_max; l++) {
    rot.M1[l] = M.block(l*l, l*l, l, l);
    rot.M2[l] = M.block(l*(l+1), l*(l+1), l+1, l+1);
  }

  return rot;
  
}
Shc::s_rot_X Shc::create_rotate_X(MatrixReal& M, int l_max) {
  
  s_rot_X rot;
  
  rot.M.resize(l_max);

  for (int l=0; l<l_max; l++) {
    rot.M[l] = M.block(l*l, l*l, 2*l+1, 2*l+1);
  }

  return rot;
  
}

// TODO: For large bands, using the full sparsity is more efficient (around l>55).
// Implementing full sparsity, the code could be optimized for high bands.
Coef Shc::compute_rotate_X2(Shc::s_rot_X2& rot, Coef& coef, int l_max) {
  
  Coef res(coef.size());
  
  // TODO HEMI-CONTIN, maxBands, etc
  VecCoef coef_pieces[4];
  for (int i=0; i<4; i++) {
    coef_pieces[i].resize(l_max);
  }

  // TODO IndexedView instead of manual copying
  for (int l=2; l<l_max; l++) {
    for (int i=0; i<4; i++) {
      coef_pieces[i][l].resize(rot.n[i][l]);
      for (int j=0; j<rot.n[i][l]; j++) {
        coef_pieces[i][l](j) = coef(l*l + rot.c[i][l](j));
      }
    }
  }
  
  res.segment(0,1) = rot.M[0][0] * coef(0);
  res.segment(1,3) = rot.M[0][1] * coef.segment(1,3);
  
  int count = 4;
  for (int l=2; l<l_max; l++) {
    for (int i=0; i<4; i++) {
      res.segment(count, rot.n[i][l]) = rot.M[i][l] * coef_pieces[i][l];
      count += rot.n[i][l];
    }
  }
  
  for (int l=2; l<l_max; l++) {
    Eigen::PermutationWrapper<VectorInt> p = PermutationWrapper<VectorInt>(rot.p[l]);
    res.segment(l*l, 2*l+1) = p * res.segment(l*l, 2*l+1).eval();
  }
  
  return res;
  
}
Shc::s_rot_X2 Shc::create_rotate_X2(MatrixReal& M, int l_max) {
  
  s_rot_X2 rot;

  for (int i=0; i<4; i++) {
    rot.c[i].resize(l_max);
    rot.r[i].resize(l_max);
    rot.n[i].resize(l_max);
    rot.M[i].resize(l_max);
  }
  
  for (int l=0; l<l_max; l++) {
    rot.n[0][l] = (l+1)/2;
    rot.n[1][l] = (l+0)/2;
    rot.n[2][l] = (l+1)/2;
    rot.n[3][l] = (l+2)/2;
  }
  
  rot.p.resize(l_max);
  for (int l=0; l<l_max; l++) {
    rot.p[l].resize(2*l+1);
  }
  
  rot.M[0][0].resize(1,1);
  rot.M[0][0] = M.block(0,0,1,1);
  
  rot.M[0][1].resize(3,3);
  rot.M[0][1] = M.block(1,1,3,3);
  
  for (int l=2; l<l_max; l++) {
    
    int ccur[4] = {0};
    int rcur[4] = {0};
    int mi, mj;
    
    for (int i=0; i<4; i++) {
      
      rot.c[i][l].resize(rot.n[i][l]);
      switch (i) {
        case 0: {mi = -l+0; break;}
        case 1: {mi = -l+1; break;}
        case 2: {mi =  l-1; break;}
        case 3: {mi =  l-0; break;}
      }
      for (mj=-l; mj<=l; mj++) {
        if (calc_rotate_isEntryZero_X(mi, mj, l) == false) {
          rot.c[i][l](ccur[i]) = mj + l;
          ccur[i]++;
        }
      }
      
      rot.r[i][l].resize(rot.n[i][l]);
      switch (i) {
        case 0: {mj =  l-1; break;}
        case 1: {mj = -l+1; break;}
        case 2: {mj = -l+0; break;}
        case 3: {mj =  l+0; break;}
      }
      for (mi=-l; mi<=l; mi++) {
        if (calc_rotate_isEntryZero_X(mi, mj, l) == false) {
          rot.r[i][l](rcur[i]) = mi + l;
          rcur[i]++;
        }
      }
      
      rot.M[i][l].resize(rot.n[i][l], rot.n[i][l]);
      for (int r=0; r<rot.n[i][l]; r++) {
        for (int c=0; c<rot.n[i][l]; c++) {
          rot.M[i][l](r,c) = M(l*l+rot.r[i][l](r), l*l+rot.c[i][l](c));
        }
      }
      
    }
    
    int count = 0;
    for (int i=0; i<4; i++) {
      for (int j=0; j<rot.n[i][l]; j++) {
        rot.p[l](count) = rot.r[i][l](j);
        count++;
      }
    }

  }
  
  return rot;
  
}

void Shc::rotate_calc_XZ(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);
    
    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      temp.segment(count,   l)   = v_rotation[index_rot_par].rotate_YP.M1[l] * in[q].coef.segment(count,   l);
      temp.segment(count+l, l+1) = v_rotation[index_rot_par].rotate_YP.M2[l] * in[q].coef.segment(count+l, l+1);
      count += 2*l+1;
    }
    
    in[q].coef = temp;
    
  }
  
}
void Shc::rotate_calc_YX(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);
    
    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      temp.segment(count, 2*l+1) = v_rotation[index_rot_par].rotate_YNXN.M[l] * in[q].coef.segment(count, 2*l+1);
      count += 2*l+1;
    }
    
    in[q].coef = temp;
    
  }
  
}
void Shc::rotate_calc_ZY(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);
    
    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      temp.segment(count, 2*l+1) = v_rotation[index_rot_par].rotate_XP.M[l] * in[q].coef.segment(count, 2*l+1);
      count += 2*l+1;
    }
    
    in[q].coef = temp;
    
  }
  
}
void Shc::rotate_calc_Z(VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par) {
  
  Coef temp;
  int count;
  
  for (int q=0; q<scp.n_total; q++) {
    if (!scp.valid[q]) {continue;}
    
    temp.resize(scp.sh_max[q][index_rot_par]);
    
    count = 0;
    for (int l=0; l<scp.l_max[q][index_rot_par]; l++) {
      if (contin_skip(in[q].contin, l)) {continue;};
      temp.segment(count, 2*l+1) = in[q].coef.segment(count, 2*l+1).array() * v_rotation[index_rot_par].rotate_Z[index].v1[l].array() + in[q].coef.segment(count, 2*l+1).reverse().array() * v_rotation[index_rot_par].rotate_Z[index].v2[l].array();
      count += 2*l+1;
    }
    
    in[q].coef = temp;
    
  }
  
}
