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

float Shc::rotation_angle(Xyz xyz) {

  MatrixRotation R = xyz2r(xyz);
  
  float result = acos((R.trace() - 1) / 2);
  
  return result;
  
}
float Shc::angular_normalization(float angle) {
  return fmod(fmod(angle, (float)(2*M_PI)) + (float)(2*M_PI), (float)(2*M_PI));
}

float Shc::hemispherical_intersection(Xyz xyz1, Xyz xyz2) {
  
  vector<Xyz> xyz;
  xyz.push_back(xyz1);
  xyz.push_back(xyz2);

  float result = hemispherical_intersection(xyz);
  
  return result;
  
}
float Shc::hemispherical_intersection(vector<Xyz> xyz) {
  
  int n = xyz.size();
  
  VectorReal rx(n);
  VectorReal ry(n);
  VectorReal rz(n);
  
  // at least two values need to be passed
  if (n <= 1)
    return 1;
  
  // calculate last row of rotation matrix
  for (int i=0; i<n; i++) {
    rx[i] = cos(xyz[i].x)*sin(xyz[i].y)*cos(xyz[i].z) + sin(xyz[i].x)*sin(xyz[i].z);
    ry[i] = cos(xyz[i].x)*sin(xyz[i].y)*sin(xyz[i].z) - sin(xyz[i].x)*cos(xyz[i].z);
    rz[i] = cos(xyz[i].x)*cos(xyz[i].y);
    if (rz[i] < 0) {
      print_warning_static("hemispherical_intersection", "rotations with negative z axis (upside-down) are not supported yet.");
      return 0;
    }
  }
  
  // calculate points of interest
  vector<float> poi;
  float num; float den;
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      num = rx[i]*rz[j]-rx[j]*rz[i];
      den = ry[i]*rz[j]-ry[j]*rz[i];
      if (den != 0) {
        poi.push_back(-atan(num/den));
        poi.push_back(-atan(num/den)+M_PI);
      } else {
        poi.push_back(M_PI/2.0);
        poi.push_back(M_PI/2.0+M_PI);
      }
    }
  }
  
  // map all points to the interval [0,2pi]
  int m = poi.size();
  for (int i=0; i<m; i++) {
    if (poi[i] < 0)
      poi[i] += 2.0*M_PI;
    if (poi[i] > 2.0*M_PI)
      poi[i] -= 2.0*M_PI;
  }
  // add points 0 and 2pi to avoid wrapping errors
  poi.push_back(0);
  poi.push_back(2.0*M_PI);
  // sort poi
  sort(poi.begin(), poi.end());
  // remove duplicate entries
  poi.erase(unique(poi.begin(), poi.end()), poi.end());
  m = poi.size();
  
  // calculate intersection
  float integral = 0;
  float q;
  float vMax; float jMax; float cur;
  for (int i=0; i<m-1; i++) {
    
    // calculate least visible hemisphere
    q = (poi[i]+poi[i+1])/2.0;
    vMax = -2.0*M_PI; jMax = 0;
    for (int j=0; j<n; j++) {
      cur = -atan( (cos(q)*rx[j]+sin(q)*ry[j]) / rz[j] );
      if (cur > vMax) {
        jMax = j; vMax = cur;
      }
    }
    
    // calculate single integrals
    int steps = (poi[i+1]-poi[i])/(2.0*M_PI)*360.0;
    VectorReal x = linspace(poi[i],poi[i+1],steps);
    VectorReal y(x.size());
    
    for (uint k=0; k<(uint)x.size(); k++) {
      float t = rx[jMax]*cos(x[k]) + ry[jMax]*sin(x[k]);
      y[k] = t / sqrt(t*t + rz[jMax]*rz[jMax]);
    }
    
    integral += trapz(x,y);
    
  }
  
  float result = (2.0*M_PI+integral) / (2.0*M_PI);

  return result;
}
float Shc::angular_difference(MatrixRotation R1, MatrixRotation R2) {
    
  MatrixRotation R = R1.transpose() * R2;
  
  // prevent floating point errors
  float t = (R.trace() - 1.0) / 2.0;
  if (t > 1) {t = 1;}
  if (t < -1) {t = -1;}
  
  // calculate angular offset
  float result = acos(t);

  return result;
  
}
float Shc::angular_difference(Xyz xyz1, Xyz xyz2) {
  
  MatrixRotation R1 = xyz2r(xyz1);
  MatrixRotation R2 = xyz2r(xyz2);
    
  return angular_difference(R1, R2);
  
}
float Shc::angular_difference(Axr axr1, Axr axr2) {
  
  MatrixRotation R1 = axr2r(axr1);
  MatrixRotation R2 = axr2r(axr2);
    
  return angular_difference(R1, R2);
  
}
float Shc::angular_difference(float angle1, float angle2) {
  
  float result;
  
  result = angular_normalization(angle1 - angle2);
  result = min(result, (float)(2*M_PI)-result);
  
  return result;
  
}
 
VectorReal Shc::linspace(float a, float b, int steps) { 
  
  if (steps < 2)
    steps = 2;
  
  VectorReal ls(steps);
  
  for (int i=0; i<steps; i++) {
    ls[i] = a + (b-a) * i/(steps-1);
  }
  
  return ls;
  
}
float Shc::trapz(VectorReal& x, VectorReal& y) {
  
  int n;
  if (x.size() == y.size()) {
    n = x.size();
  } else {
    print_warning_static("trapz", "unequal vector lengths");
    return 0;
  }
  
  float sum = 0;
  for (int i=0; i<n-1; i++) {
    sum += (x[i+1]-x[i])*(y[i+1]+y[i]);
  }
  
  sum /= 2;
  
  return sum;
  
}

bool Shc::tick() {
  
  if (bool_tick == false) {
    clock_gettime(CLOCK_MONOTONIC, &measureTick);
    bool_tick = true;
    return true;
  }
  
  print_warning("tick", "A previously started measurement is not finished yet.");
  return false;
  
}
float Shc::tock() {
  
  if (bool_tick == true) {
  
    struct timespec measureTock;
    clock_gettime(CLOCK_MONOTONIC, &measureTock);
    float ms = (measureTock.tv_sec * 1.0e3 + measureTock.tv_nsec / 1.0e6) - (measureTick.tv_sec * 1.0e3 + measureTick.tv_nsec / 1.0e6);
    bool_tick = false;
    return ms;
    
  }
  
  print_warning("tick", "No measurement has been started yet.");
  return -1;
  
}

VectorReal Shc::get_feature(Shpm& shpm, e_feature feature) {
  
  if (!check_init() || !check_init(shpm)) {return VectorReal();}
  
  VectorReal result;
  
  switch (feature) {
    case ISE:
      result = calc_feature_ISE(shpm);
      break;
    case AS:
      result = calc_feature_amplitudespectrum(shpm);
      break;
    case BS:
      result = calc_feature_bispectrum(shpm);
      break;
    case BS_DENSE:
      result = calc_feature_bispectrum_dense(shpm);
      break;
  }
  
  return result;
  
}
VectorReal Shc::get_feature_difference(VecShpm& vshpm1, Shpm& shpm2, e_feature feature) {
  
  int n = vshpm1.size();
  VectorReal result(n);
  
  if (!check_init() || !check_init(vshpm1) || !check_init(shpm2)) {
    result.setZero();
    return result;
  }
  
  VecShpm s1(n);
  VecShpm s2(n);
  
  for (int i=0; i<n; i++) {
    if (trifold(vshpm1[i].contin) != trifold(shpm2.contin)) {
      s1[i] = contin_convert(vshpm1[i], shpm2.contin);
    } else {
      s1[i] = vshpm1[i];
    }
  }
  
  switch (tangent_distance_mode) {
    case NONE:
      for (int i=0; i<n; i++) {
        s2[i] = shpm2;
      }
      break;
    case ONESIDED:
      s2 = tangent_distance_calc(vshpm1, shpm2);
      break;
    case TWOSIDED:
      for (int i=0; i<n; i++) {
        tangent_distance_calc(s1[i], shpm2, s1[i], s2[i]);
      }
      break;
  }
  
  for (int i=0; i<n; i++) {
    switch (feature) {
      case ISE:
        result[i] = calc_difference_ISE(s1[i], s2[i]);
        break;
      case AS:
        result[i] = calc_difference_AS(s1[i], s2[i]);
        break;
      case BS:
        result[i] = calc_difference_BS(s1[i], s2[i]);
        break;
      case BS_DENSE:
        result[i] = calc_difference_BS_dense(s1[i], s2[i]);
        break;
    }
  }
  
  return result;
  
}
float Shc::get_feature_difference(Shpm& shpm1, Shpm& shpm2, e_feature feature) {
  
  if (!check_init()) {return 0;}
  
  VecShpm vshpm1;
  vshpm1.push_back(shpm1);
  
  VectorReal result = get_feature_difference(vshpm1, shpm2, feature);
  
  return result[0];
  
}

int Shc::get_feature_index_sh(Shpm& shpm, int l, int m) {

  if (!check_init()) {return 0;}
  
  if (l < 0 || l > shpm.l_max || abs(m) > abs(l)) {
    print_warning("get_index_sh", "invalid indices / indices out of bounds");
    return -1;
  }
  
  int i = l2sh(l, shpm.contin);
  i += l + m;
  
  return i;
   
}
int Shc::get_feature_index_as(Shpm& shpm, int l) {

  if (!check_init()) {return 0;}
  
  if (l < 0 || l > shpm.l_max) {
    print_warning("get_index_sh", "invalid indices / indices out of bounds");
    return -1;
  }
  
  int countRes = 0;
  for (int il=0; il<shpm.l_max; il++) {
    
    if (contin_skip(shpm.contin, il) == true) {continue;}

    if (il == l) {
      return countRes;
    }
    
    countRes++;
    
  } 
  
  return -1;
   
}
int Shc::get_feature_index_bs(Shpm& shpm, int l1, int l2, int i) {
  
  if (!check_init()) {return 0;}
  
  if (l1 < 0 || l1 > shpm.l_max || l2 < 0 || l2 > shpm.l_max || i < 0 || i > shpm.l_max) {
    print_warning("get_index_sh", "invalid indices / indices out of bounds");
    return -1;
  }
  
  int l_max = min(n_bands_CG, shpm.l_max); 
  
  int countRes = 0;
  for (int il1=0; il1<l_max; il1++) {
    for (int il2=il1; il2<l_max; il2++) {
      
      if (contin_skip(shpm.contin, il1) || contin_skip(shpm.contin, il2)) {continue;}
      
      for (int ii=il2-il1; ii<=il2+il1; ii++) {
        
        if (contin_skip(shpm.contin, ii)) {continue;}

        if (il1 == l1 && il2 == l2 && ii == i) {
          return countRes;
        }
          
        countRes++;
        
      }

    }
  }
  
  return -1;
   
}
int Shc::get_feature_index_bs_dense(Shpm& shpm, int l1, int l2, int i) {
  
  if (!check_init()) {return 0;}
  
  if (l1 < 0 || l1 > shpm.l_max || l2 < 0 || l2 > shpm.l_max || i < 0 || i > shpm.l_max) {
    print_warning("get_index_sh", "invalid indices / indices out of bounds");
    return -1;
  }
  
  if (trifold(shpm.contin) != FULL) {
    print_warning("calc_feature_bispectrum_dense", "Only Shpm without hemispherical continuation (i.e. contin == FULL) can be used.");
    return -1;
  }
  
  int l_max = min(n_bands_CG, shpm.l_max); 
  
  int countRes = 0; int il2; int ii;
  for (int il1=1; il1<l_max; il1++) {
      
    il2 = il1-1;
    ii = 1;
    
    if (il1 == l1 && il2 == l2 && ii == i) {
      return countRes;
    }
    
    countRes++;

  }
  
  return -1;
   
}



void Shc::get_surface_data(VectorReal& theta, VectorReal& phi, VectorReal& weight) {
  
  if (!check_init()) {return;}
  
  int n = lut_surf.n_angles;
  
  theta.resize(n);
  phi.resize(n);
  weight.resize(n);
  
  for (int i=0; i<n; i++) {
    theta[i]  = lut_surf.unit[i].sphericalCoor.theta;
    phi[i]    = lut_surf.unit[i].sphericalCoor.phi;
    weight[i] = lut_surf.unit[i].weight;
  }
  
}
int Shc::get_surface_size() {
  
  if (!check_init()) {return 0;}
  
  return lut_surf.n_angles;
  
}
int Shc::get_surface_nearest(float theta, float phi) {
  
  if (!check_init()) {return 0;}
  
  if (bool_sphere_mode == false) {
    print_warning("get_surface_nearest", "Closest surface points can only be calculated using non-custom sample points");
    return 0;
  }
  
  s_sphericalCoor sphericalCoor;
  sphericalCoor.theta = theta;
  sphericalCoor.phi   = phi;
  
  vector<int> index = get_nearest_indices(sphericalCoor);
  
  float wt, wp, d;
  float min_d = 1e12;
  int min_index = 0;
  int n = index.size();
  
  for (int i=0; i<n; i++) {
    wt = angular_difference(sphericalCoor.theta, lut_surf.unit[index[i]].sphericalCoor.theta);
    wp = angular_difference(sphericalCoor.phi,   lut_surf.unit[index[i]].sphericalCoor.phi);
    d = sqrt(wt*wt+wp*wp);
    if (d < min_d) {
      min_d = d;
      min_index = i;
    }
  }
  
  return index[min_index];
  
}
int Shc::get_surface_nearest(Coor3d coor) {
  
  if (!check_init()) {return 0;}
  
  s_coor3 coor3(coor);
  s_sphericalCoor sphericalCoor = coor32sphericalCoor(coor3);

  return get_surface_nearest(sphericalCoor.theta, sphericalCoor.phi);
  
}

Xyz Shc::transpose(Xyz xyz) {
  
  MatrixRotation R = xyz2r(xyz);
  R.transposeInPlace();
  
  Xyz result = r2xyz(R);
  
  return result;
  
}
Axr Shc::transpose(Axr axr) {
  
  Axr result;
  
  result.axis  = axr.axis;
  result.angle = -axr.angle;
  
  return result;
  
}
MatrixRotation Shc::transpose(MatrixRotation r) {
  
  MatrixRotation R = r.transpose();
  
  return R;
  
}

Xyz Shc::product(Xyz xyz1, Xyz xyz2) {

  MatrixRotation R1 = xyz2r(xyz1);
  MatrixRotation R2 = xyz2r(xyz2);
  MatrixRotation R  = R1 * R2;
  
  Xyz result = r2xyz(R);
  
  return result;
  
}
Axr Shc::product(Axr axr1, Axr axr2) {

  MatrixRotation R1 = axr2r(axr1);
  MatrixRotation R2 = axr2r(axr2);
  MatrixRotation R  = R1 * R2;
  
  Axr result = r2axr(R);
  
  return result;
  
}
MatrixRotation Shc::product(MatrixRotation r1, MatrixRotation r2) {

  MatrixRotation R = r1 * r2;
  
  return R;
  
}

void Shc::get_bands(int& n_bands, int& n_bands_CG) {
  
  if (!check_init()) {
    n_bands    = -1;
    n_bands_CG = -1;
  } else {
    n_bands    = this->n_bands;
    n_bands_CG = this->n_bands_CG;
  }
  
}
