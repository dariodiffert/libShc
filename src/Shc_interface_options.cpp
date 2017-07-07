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

void Shc::set_feature_norm(e_feature_norm feature_norm) {
  
  if (!check_init()) {return;}
  
  this->feature_norm = feature_norm;
  
}
e_feature_norm Shc::get_feature_norm() {
  
  if (!check_init()) {}
  
  return this->feature_norm;
  
}
void Shc::set_feature_normalize(bool feature_normalize) {
  
  if (!check_init()) {return;}
  
  this->bool_feature_normalize = feature_normalize;
  
}
bool Shc::get_feature_normalize() {
  
  if (!check_init()) {}
  
  return this->bool_feature_normalize;
  
}
void Shc::set_print_level(e_print_level print_level) {
  this->print_level = print_level;
}
e_print_level Shc::get_print_level() {
  return print_level;
}
void Shc::set_filter(e_filter filter) {
  
  if (!check_init()) {return;}
  
  switch (filter) {
    case NEAREST:
      callback_filter_mat = &Shc::filter_nearest_mat;
      callback_filter_raw = &Shc::filter_nearest_raw;
      break;
    case BILINEAR:
      callback_filter_mat = &Shc::filter_bilinear_mat;
      callback_filter_raw = &Shc::filter_bilinear_raw;
      break;
  }
  
  this->filter = filter;
}
e_filter Shc::get_filter() {
  
  if (!check_init()) {}
  
  return filter;
  
}
void Shc::set_linearize_translations(bool linearize) {
  
  if (!check_init()) {return;}
  
  bool_linearize_translations = linearize;
  
}
bool Shc::get_linearize_translations() {
  
  if (!check_init()) {}
  
  return bool_linearize_translations;
  
}

void Shc::set_noise_mask(float angle) {
  
  if (!check_init()) {return;}
  
  if (angle < 0 || angle > 2.0*M_PI + 1e-6) {
    print_warning("set_noise_mask", "The angle has to be in the interval [0, 2*PI]");
  }
  
  Surf surf(lut_surf.n_angles);
  
  for (int i=0; i<lut_surf.n_angles; i++) {
    if (lut_surf.unit[i].sphericalCoor.theta > angle / 2.0) {
      surf[i] = 1;
    } else {
      surf[i] = 0;
    }
  }
  
  set_noise_mask(surf);
  
}
void Shc::set_noise_mask(MatrixReal& mask) {
  
  if (!check_init()) {return;}
  
  Surf surf = matrix2surf(mask);
  set_noise_mask(surf);
  
}
void Shc::set_noise_mask(Surf& surf) {
  
  if (!check_init()) {return;}
  
  int n = surf.size();
  
  if (n != lut_surf.n_angles) {
    print_warning("set_noise_mask","size of surf does not match the surface parameters of the Shc instance");
    return;
  }
  
  for (int i=0; i<n; i++) {
    if (surf(i) < 0) {surf(i) = 0;}
    if (surf(i) > 1) {surf(i) = 1;}
  }

  noise_mask     = surf;
  noise_mask_sum = surf.sum();
  
}
Surf Shc::get_noise_mask() {
  
  if (!check_init()) {return Surf();}
  
  return noise_mask;
  
}

void Shc::set_noise_enable(bool enable) {
  
  if (!check_init()) {return;}
  
  noise_enable = enable;
  
}
bool Shc::get_noise_enable() {
  
  if (!check_init()) {}
  
  return noise_enable;
  
}

void Shc::set_noise_amplifier(float amplifier) {
  
  if (!check_init()) {return;}
  
  noise_amplifier = amplifier;
  
}
float Shc::get_noise_amplifier() {
  
  if (!check_init()) {}
  
  return noise_amplifier;
  
}





