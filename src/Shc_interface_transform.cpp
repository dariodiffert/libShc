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
 
MatrixReal Shc::create_matrix_rotation(Xyz xyz, int l_max) {
 
  if (!check_init()) {return MatrixReal();}
   
  MatrixReal X = calc_rotate_axis(xyz.x, AXIS_X, l_max);
  MatrixReal Y = calc_rotate_axis(xyz.y, AXIS_Y, l_max);
  MatrixReal Z = calc_rotate_axis(xyz.z, AXIS_Z, l_max);
  
  MatrixReal result = X*Y*Z;
  
  return result;
  
}
MatrixReal Shc::create_matrix_warp(Coor3d translation, e_trans_type trans_type, e_contin contin, int l_max) {
  
  if (!check_init()) {return MatrixReal();}

  s_coor3 trans(translation);
  
  Xyz rot;
  
  if (l_max < 0 || l_max > n_bands) {
    l_max = n_bands;
  }

  MatrixReal M = calc_translation_simple(translation, trans_type, contin, l_max); 
  
  return M;
  
}

void Shc::set_tangent_distance(e_tangent_distance_mode mode, float spring_constant) {
  
  if (!check_init()) {return;}
  
  this->tangent_distance_mode            = mode;
  this->tangent_distance_spring_constant = spring_constant;
  
}
void Shc::get_tangent_distance(e_tangent_distance_mode& mode, float& spring_constant) {
  
  if (!check_init()) {return;}
  
  mode            = this->tangent_distance_mode;
  spring_constant = this->tangent_distance_spring_constant;
  
}
int Shc::add_tangent_distance(MatrixReal& matrix_transform, e_matrix_type matrix_type, float tolerance, int l_max) {
  
  if (!check_init()) {return false;}
  
  int res = intern_add_transform(tangent_distance, matrix_transform, matrix_type, tolerance, l_max);
  
  return res;
    
}
int Shc::add_tangent_distance_rotation(e_axis axis, float angle, int l_max) {
  
  if (!check_init()) {return 0;}
  
  if (l_max > n_bands) {
    print_warning("add_tangent_distance", "l_max is chosen higher than the maximal number of bands.");
    l_max = n_bands;
  }
  
  MatrixReal matP = calc_rotate_axis( angle, axis, n_bands);
  MatrixReal matN = calc_rotate_axis(-angle, axis, n_bands);
  MatrixReal matE = MatrixXf::Identity(n_bands*n_bands,n_bands*n_bands);
  MatrixReal mat = matP - matN + matE;
  
  int result;
  result = add_tangent_distance(mat, BANDWISE, 0, l_max);
  
  return result;
  
}
int Shc::add_tangent_distance_translation(e_axis axis, float dist, e_contin contin, int l_max) {
  
  if (!check_init()) {return 0;}
  
  if (l_max > n_bands) {
    print_warning("add_tangent_distance", "l_max is chosen higher than the maximal number of bands.");
    l_max = n_bands;
  }
  
  MatrixReal matP, matN;
  
  switch (axis) {
    case AXIS_X:
      matP = create_matrix_warp(Coor3d( dist,0,0), VISUAL, contin);
      matN = create_matrix_warp(Coor3d(-dist,0,0), VISUAL, contin);
      break;
    case AXIS_Y:
      matP = create_matrix_warp(Coor3d(0, dist,0), VISUAL, contin);
      matN = create_matrix_warp(Coor3d(0,-dist,0), VISUAL, contin);
      break;
    case AXIS_Z:
      matP = create_matrix_warp(Coor3d(0,0, dist), VISUAL, contin);
      matN = create_matrix_warp(Coor3d(0,0,-dist), VISUAL, contin);
      break;
  }
  
  MatrixReal matE = MatrixXf::Identity(n_bands*n_bands,n_bands*n_bands);
  MatrixReal mat = matP - matN + matE;
  
  int result;
  result = add_tangent_distance(mat, SPARSE, 1e-6, l_max);
  
  return result;
  
}
void Shc::clear_tangent_distance() {
  
  if (!check_init()) {return;}
  
  intern_clear_transform(tangent_distance);

}
void Shc::compass_configure_tangent_distance(int index, VectorInt td_active) {
  // index: Index of the compass step (each call to init_rotations_sphere/cone creates one compass step)
  // td_active: Contains the indices of tangent distance transformations applied for the compass step. If the vector is empty, all tangent distance transforms are applied.
  
  if (!check_init()) {return;}
  
  if ((index < 0) || (index >= n_rot_par)) {
    stringstream ss;
    ss << "The specified index (" << index << ") does not correspond to an initialized rotation (rotations indices 0, ...," << n_rot_par << ")";
    print_warning("compass_configure_tangent_distance", ss.str());
    return;
  }
  
  v_rotation[index].td_active = td_active;
  
}
int Shc::get_tangent_distance_size() {
  
  if (!check_init()) {return 0;}
  
  return tangent_distance.unit.size();
}

int Shc::get_transform_size() {
  
  if (!check_init()) {return 0;}
  
  return transforms.unit.size();
}
int Shc::add_transform(MatrixReal& matrix_transform, e_matrix_type matrix_type, float tolerance, int l_max) {
  
  if (!check_init()) {return 0;}
  
  return intern_add_transform(transforms, matrix_transform, matrix_type, tolerance, l_max);
    
}
void Shc::clear_transform() {
  
  if (!check_init()) {return;}
  
  intern_clear_transform(transforms);
  
}

bool Shc::save_tangent_distance(string file) {
  
  stringstream ss;
  ss << "save_tangent_distance: '" << file << "'\n";
  print(ss.str());
  
  if (!check_init()) {return false;}
  
  bool result = intern_save_transform(tangent_distance, file);
  if (result == false) {  
    print_warning("save_tangent_distance","Could not create file.");
  }
  
  return result;
  
}
bool Shc::load_tangent_distance(string file) {
  
  stringstream ss;
  ss << "load_tangent_distance: '" << file << "'\n";
  print(ss.str());
  
  if (!check_init()) {return false;}
  
  clear_tangent_distance();
  bool result = intern_load_transform(tangent_distance, file);
  if (result == false) {  
    print_warning("load_tangent_distance","Could not create file.");
  }
  
  return result;
  
}
bool Shc::save_transform(string file) {
  
  stringstream ss;
  ss << "save_transform: '" << file << "'\n";
  print(ss.str());
  
  if (!check_init()) {return false;}
  
  bool result = intern_save_transform(transforms, file);
  if (result == false) {  
    print_warning("save_transform","Could not create file.");
  }
  
  return result;
  
}
bool Shc::load_transform(string file) {
  
  stringstream ss;
  ss << "load_transform: '" << file << "'\n";
  print(ss.str());
  
  if (!check_init()) {return false;}
  
  clear_transform();
  bool result = intern_load_transform(transforms, file);
  if (result == false) {  
    print_warning("load_transform","Could not load file.");
  }
  
  return result;
  
}


