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

void Shc::create_fast_indices(VectorReal angles, float resolution, s_rotation &rotation) {
  
  for (uint j=0; j<(uint)rotation.psi.size(); j++) {
    rotation.fast_indices_z.push_back(j);
  }
  
  for (uint i=0; i<(uint)angles.size(); i++) {
    for (uint j=0; j<(uint)rotation.psi.size(); j++) {
      if (angular_difference(rotation.psi[j], angles[i]) < resolution / 2) {
        rotation.fast_indices_y.push_back(j);
        rotation.fast_indices_x.push_back(j);
        break;
      }
    }
  }
  rotation.n_fast_indices_z = rotation.fast_indices_z.size();
  rotation.n_fast_indices_y = rotation.fast_indices_y.size();
  rotation.n_fast_indices_x = rotation.fast_indices_x.size();
  
}

void Shc::index2sh(int i, int& l, int& m) {
  l = lut_index2sh[i].l;
  m = lut_index2sh[i].m;
}
int Shc::sh2index(int l, int m) {
  return ((l+1)*l+m);
}
int Shc::sha2index(s_surf& surface, int l, int m, int a) {
  return sh2index(l,m)*surface.n_angles+a;
}
int Shc::l2sh(int l, e_contin contin) {
  
  switch (contin) {
    case HEMI_RM:
      l = ceil(l/2.0);
      return l*(2*l-1);
    case HEMI_RMN:
      l = ceil((l+1)/2.0);
      return l*(2*l-3)+2;
    default:
      return l*l;
  }
  
}

bool Shc::contin_skip(e_contin contin, int l) {
  
  if (l == 0)
    return false;
  
  if (contin == HEMI_RM  && l%2 == 1)
    return true;
  
  if (contin == HEMI_RMN && l%2 == 0)
    return true;
  
  return false;
  
}
bool Shc::even(int n) {
  if (abs(n)%2 == 0) {
    return true;
  } else {
    return false;
  }
}
bool Shc::odd(int n) {
  if (abs(n)%2 == 1) {
    return true;
  } else {
    return false;
  }
}

e_contin Shc::trifold(e_contin contin) {
  switch (contin) {
    case HEMI_RM:
      return HEMI_RM;
    case HEMI_RMN:
      return HEMI_RMN;
    default:
      return FULL;
  }
}
