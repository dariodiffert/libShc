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
 
Xyz Shc::r2xyz(MatrixRotation r) {
  
  Xyz xyz;

  if (r(0,2) < 1) {
    if (r(0,2) > -1) {
      xyz.y = asin(r(0,2));
      xyz.x = atan2(-r(1,2),r(2,2));
      xyz.z = atan2(-r(0,1),r(0,0));
    } else {
      xyz.y = -M_PI/2;
      xyz.x = -atan2(-r(1,0),r(1,1));
      xyz.z = 0;
    }
  } else {
    xyz.y = M_PI/2;
    xyz.x = atan2(r(1,0),r(1,1));
    xyz.z = 0;
  }
  
  return xyz;
  
}
MatrixRotation Shc::xyz2r(Xyz xyz) {
  
  float x = xyz.x;
  float y = xyz.y;
  float z = xyz.z;
  
  MatrixRotation R;
  
  R(0,0) = cos(z)*cos(y);  
  R(1,0) = sin(z)*cos(x)+cos(z)*sin(y)*sin(x);       
  R(2,0) = sin(z)*sin(x)-cos(z)*sin(y)*cos(x); 
  
  R(0,1) = -sin(z)*cos(y); 
  R(1,1) =  cos(z)*cos(x)-sin(z)*sin(y)*sin(x); 
  R(2,1) =  cos(z)*sin(x)+sin(z)*sin(y)*cos(x); 
  
  R(0,2) =  sin(y);
  R(1,2) = -cos(y)*sin(x);
  R(2,2) =  cos(y)*cos(x);
  
  return R;
  
}
Xyz Shc::axr2xyz(Axr axr) {
  
  MatrixRotation R = axr2r(axr);
  Xyz xyz = r2xyz(R);
  
  return xyz;
  
}
MatrixRotation Shc::axr2r(Axr axr) {
  
  Coor3d rot_axis = axr.axis / axr.axis.norm();
  
  MatrixRotation K;
  K <<            0, -rot_axis(2),  rot_axis(1),
        rot_axis(2),            0, -rot_axis(0),
       -rot_axis(1),  rot_axis(0),            0;
       
  MatrixRotation I = MatrixRotation::Identity(3, 3);
  MatrixRotation R = I + sin(axr.angle) * K + (1-cos(axr.angle)) * K * K;
  
  return R;
  
}
Axr Shc::r2axr(MatrixRotation r) {
  
  MatrixRotation T = r - r.transpose();
  
  Axr axr;
  axr.axis(0) = -T(1,2);
  axr.axis(1) =  T(0,2);
  axr.axis(2) = -T(0,1);
  
  float n = axr.axis.norm();
  
  if (abs(n) > 1e-8) {
    axr.axis /= n;
  } else {
    axr.axis = Coor3d(1,0,0);
  }
  
  axr.angle = acos((r.trace() - 1) / 2);
  
  return axr;
  
}
Axr Shc::xyz2axr(Xyz xyz) {

  MatrixRotation R = xyz2r(xyz);
  Axr axr = r2axr(R);
  
  return axr;
  
}

MatrixRotation Shc::tilt2r(float tilt_dir, float tilt_angle) {
  // creates a rotation which tilts the z axis by tilt_angle degree into the direction tilt_dir in [0,2pi]
  
  MatrixRotation RZ = xyz2r(Xyz(0,0,tilt_dir));
  
  Coor3d rot_axis = RZ * Coor3d(0,1,0);
  
  Axr axr(tilt_angle, rot_axis);
  MatrixRotation RT = axr2r(axr);
  
  return RT;
  
}
Xyz Shc::tilt2xyz(float tilt_dir, float tilt_angle) {

  MatrixRotation RT = tilt2r(tilt_dir, tilt_angle);
  Xyz xyz = r2xyz(RT);
  
  return xyz;
  
}
Axr Shc::tilt2axr(float tilt_dir, float tilt_angle) {
  
  MatrixRotation RT = tilt2r(tilt_dir, tilt_angle);
  Axr axr = r2axr(RT);
  
  return axr;
  
}

void Shc::vector2sphericalCoordinate(Coor3d vector, float& theta, float& phi) {
  
  s_coor3 coor3(vector);
  s_sphericalCoor sc = coor32sphericalCoor(coor3);
  
  theta = sc.theta;
  phi   = sc.phi;
  
}
Coor3d Shc::sphericalCoordinate2vector(float theta, float phi) {
  
  s_sphericalCoor sc;
  sc.theta = theta;
  sc.phi   = phi;
  
  s_coor3 coor3 = sphericalCoor2coor3(sc);
  
  Coor3d vector;
  vector(0) = coor3.x;
  vector(1) = coor3.y;
  vector(2) = coor3.z;
  
  return vector;
  
}

Surf Shc::matrix2surf(MatrixReal& matrix) {
  
  if (!check_init()) {return Surf();}
  
  int h = matrix.rows();
  int w = matrix.cols();
  float x; float y;
  Surf surf(lut_surf.n_angles);  
  
  if (w==0 || h==0) {
    print_warning("matrix2surf", "empty matrix passed as argument");
    surf.setZero();
    return surf;
  }
  
  for (int i=0; i<lut_surf.n_angles; i++) {
    y =         lut_surf.unit[i].sphericalCoor.theta / M_PI     * (h-1);
    x = (w-1) - lut_surf.unit[i].sphericalCoor.phi   / (2*M_PI) * (w-1);
    surf[i] = callback_filter_mat(matrix, h, w, y, x);
  }
  
  return surf;
  
}
Surf Shc::matrix2surf(unsigned char* data, int width, int height) {
  
  if (!check_init()) {return Surf();}
  
  float x; float y;
  Surf surf(lut_surf.n_angles);  
  
  for (int i=0; i<lut_surf.n_angles; i++) {
    y =             lut_surf.unit[i].sphericalCoor.theta / M_PI     * (height-1);
    x = (width-1) - lut_surf.unit[i].sphericalCoor.phi   / (2*M_PI) * (width-1);
    surf[i] = callback_filter_raw(data, height, width, y, x);
  }
  
  return surf;
  
}
Shpm Shc::matrix2shpm(MatrixReal& matrix, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  if (matrix.rows() == 0 || matrix.cols() == 0) {
    print_warning("matrix2shpm", "empty matrix passed as argument");
    return Shpm();
  }
  
  if (l_max > n_bands) {
    l_max = n_bands;
    print_warning("matrix2shpm", "l_max set to a higher band than n_bands. l_max will be adjusted");
  }
  if (l_max == -1) {
    l_max = n_bands;
  }
  
  Surf surf = matrix2surf(matrix);
  Shpm shpm = surf2shpm(surf, contin, l_max);
  
  return shpm;
  
}
Shpm Shc::matrix2shpm(unsigned char* data, int width, int height, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  
  if (l_max > n_bands) {
    l_max = n_bands;
    print_warning("matrix2shpm", "l_max set to a higher band than n_bands. l_max will be adjusted");
  }
  if (l_max == -1) {
    l_max = n_bands;
  }
  
  Surf surf = matrix2surf(data, width, height);
  Shpm shpm = surf2shpm(surf, contin, l_max);
  
  return shpm;
  
}
MatrixReal Shc::shpm2matrix(Shpm& shpm) {
  
  if (!check_init() || !check_init(shpm)) {return MatrixReal();}

  MatrixReal out;
  
  Surf surf;  
  surf = intern_shpm2output(shpm);
  
  out = intern_output2matrix(surf);
  
  return out;
  
}

Shpm Shc::shpm2shpm(Shpm& shpm, e_contin contin, int l_max) {
  
  if (!check_init() || !check_init(shpm)) {return Shpm();}
  
  Shpm result = contin_convert(shpm, contin, l_max);
  
  return result;
  
}
Shpm Shc::coef2shpm(Coef& coef, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  
  if (l_max < 0) {
    l_max = n_bands;
  }
  if (l_max > n_bands) {
    print_warning("create_shpm", "l_max parameter exceeds the maximal number of bands (n_bands) which can be used by this instance of Shc and is adjusted.");
    l_max = n_bands;
  }
  
  
  int m = coef.size();
  int n = min(m, l2sh(l_max, FULL));
  
  Shpm shpm;
  shpm.coef.resize(l2sh(l_max, FULL));
  shpm.coef.setZero();  
  shpm.coef.segment(0,n) = coef.segment(0,n);
  shpm.contin = FULL;
  shpm.l_max  = l_max;
  
  Shpm result = contin_convert(shpm, contin, l_max);
  
  return shpm;
  
}
Shpm Shc::surf2shpm(Surf& surf, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  
  if (l_max < 0) {
    l_max = n_bands;
  }
  if (l_max > n_bands) {
    print_warning("create_shpm", "l_max parameter exceeds the number of n_bands and is adjusted.");
    l_max = n_bands;
  }
  
  Shpm result;
  stringstream ss;
  
  if (surf.size() == lut_surf.n_angles) {
    result = intern_surf2shpm(surf, contin, l_max, true, contin);
  } else {
    ss << "input size (" << surf.size() << ") ";
    ss << "does not match the number of surface points (" << lut_surf.n_angles << ").";
    print_warning("create_shpm", ss.str());
  }
  
  return result;
  
}
Coef Shc::shpm2coef(Shpm& shpm) {
  
  if (!check_init() || !check_init(shpm)) {return VectorReal();}
  
  Shpm full = contin_convert(shpm, FULL);
  
  return full.coef;

}
Surf Shc::shpm2surf(Shpm& shpm) {
  
  if (!check_init() || !check_init(shpm)) {return VectorReal();}
  
  return intern_shpm2surf(shpm);
  
}

MatrixReal Shc::ocamcalib2matrix(MatrixReal& matrix) {

  if (!check_init()) {return MatrixReal();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2matrix", "OCamCalib is not initialized!");
    return MatrixReal();
  }
  
  if (ocam_model.height != matrix.rows() || ocam_model.width != matrix.cols()) {
    print_warning("ocamcalib2matrix", "Matrix does not have the dimensions specified by the loaded ocamcalib calibration file!");
    return MatrixReal();
  }
  
  MatrixReal mat(output_height, output_width);
  mat.setZero();
  
  s_coor2 coor2;
  for (float y=0; y<output_height; y++) {
    for (float x=0; x<output_width; x++) {
      mat(y,x) = callback_filter_mat(matrix, ocam_model.height, ocam_model.width, ocam_model.mappingConversionY(y,x), ocam_model.mappingConversionX(y,x));
    }
  }
  
  return mat;
  
}
MatrixReal Shc::ocamcalib2matrix(unsigned char* data) {

  if (!check_init()) {return MatrixReal();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2matrix", "OCamCalib is not initialized!");
    return MatrixReal();
  }
  
  MatrixReal mat(output_height, output_width);
  mat.setZero();
  
  s_coor2 coor2;
  for (float y=0; y<output_height; y++) {
    for (float x=0; x<output_width; x++) {
      mat(y,x) = callback_filter_raw(data, ocam_model.height, ocam_model.width, ocam_model.mappingConversionY(y,x), ocam_model.mappingConversionX(y,x));
    }
  }
  
  return mat;
  
}
void Shc::ocamcalib2matrix(unsigned char* data, unsigned char* data_out) {

  if (!check_init()) {return;}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2matrix", "OCamCalib is not initialized!");
    return;
  }
  
  int i=0;
  s_coor2 coor2;
  float res;
  for (float y=0; y<output_height; y++) {
    for (float x=0; x<output_width; x++) {
      res = callback_filter_raw(data, ocam_model.height, ocam_model.width, ocam_model.mappingConversionY(y,x), ocam_model.mappingConversionX(y,x)) * 255.0;
      if (res > 255) {res = 255;}
      if (res < 0) {res = 0;}
      data_out[i] = round(res);
      i++;
    }
  }
  
}
Surf Shc::ocamcalib2surf(MatrixReal& matrix) {
  
  if (!check_init()) {return Surf();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2surf", "OCamCalib is not initialized!");
    return Surf();
  }
  
  if (ocam_model.height != matrix.rows() || ocam_model.width != matrix.cols()) {
    print_warning("ocamcalib2surf", "Matrix does not have the dimensions specified by the loaded ocamcalib calibration file!");
    return Surf();
  }
  
  VectorReal surf(lut_surf.n_angles);  
  
  for (int i=0; i<lut_surf.n_angles; i++) {
    surf[i] = callback_filter_mat(matrix, ocam_model.height, ocam_model.width, ocam_model.mapped[i].y, ocam_model.mapped[i].x);
  }
  
  return surf;
  
}
Surf Shc::ocamcalib2surf(unsigned char* data) {
  
  if (!check_init()) {return Surf();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2surf", "OCamCalib is not initialized!");
    return Surf();
  }
  
  VectorReal surf(lut_surf.n_angles);  
  
  for (int i=0; i<lut_surf.n_angles; i++) {
    surf[i] = callback_filter_raw(data, ocam_model.height, ocam_model.width, ocam_model.mapped[i].y, ocam_model.mapped[i].x);
  }
  
  return surf;
  
}
Shpm Shc::ocamcalib2shpm(MatrixReal& matrix, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2shpm", "OCamCalib is not initialized!");
    return Shpm();
  }
  
  Surf surf = ocamcalib2surf(matrix);

  Shpm shpm = surf2shpm(surf, contin, l_max);
  
  return shpm;
  
}
Shpm Shc::ocamcalib2shpm(unsigned char* data, e_contin contin, int l_max) {
  
  if (!check_init()) {return Shpm();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib2shpm", "OCamCalib is not initialized!");
    return Shpm();
  }
  
  Surf surf = ocamcalib2surf(data);

  Shpm shpm = surf2shpm(surf, contin, l_max);
  
  return shpm;
  
}

Coor3d Shc::ocamcalib_cam2world(Coor2d coor2) {
  
  if (!check_init()) {return Coor3d();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib_cam2world", "OCamCalib is not initialized!");
    return Coor3d();
  }
  
  s_coor2 c2(coor2);
  s_coor3 c3;
  
  cam2world(c2, c3);
  
  Coor3d result;
  result(0) = c3.x;
  result(1) = c3.y;
  result(2) = c3.z;
  
  return result;
  
}
Coor2d Shc::ocamcalib_world2cam(Coor3d coor3) {
  
  if (!check_init()) {return Coor2d();}
  
  if (bool_init_ocamcalib == false) {
    print_warning("ocamcalib_world2cam", "OCamCalib is not initialized!");
    return Coor2d();
  }
  
  s_coor3 c3(coor3);
  s_coor2 c2;
  
  world2cam(c2, c3);
  
  Coor2d result;
  result(0) = c2.x;
  result(1) = c2.y;
  
  return result;
  
}

VecShpm Shc::pointcloud2shpm(Pointcloud& pointcloud, e_contin contin, int l_max) {

  if (!check_init()) {return VecShpm();}
  
  // create surfaces, where the point cloud is projected on
  int n_surface = get_surface_size();
  vector<Surf> v_surf(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    v_surf[i].resize(n_surface);
    v_surf[i].setZero();
  }
  
  // project pointcloud on slices (linearly interpolated)
  int n_pointcloud = pointcloud.size();
  for (int k=0; k<n_pointcloud; k++) {
    
    // find closest point on the surface
    int j = get_surface_nearest(pointcloud[k]);
    
    if (pointcloud[k].norm() <= lut_translate.slices[0]) {
      v_surf[0](j)++;
    } else if (pointcloud[k].norm() >= lut_translate.slices[lut_translate.n_slices-1]) {
      v_surf[lut_translate.n_slices-1](j)++;
    } else {
      for (int i=0; i<lut_translate.n_slices-1; i++) {
        
        if (pointcloud[k].norm() >= lut_translate.slices(i) && pointcloud[k].norm() < lut_translate.slices(i+1)) {
          float d1 = pointcloud[k].norm() - lut_translate.slices(i);
          float d2 = lut_translate.slices(i+1) - pointcloud[k].norm();
          float d = d1+d2;
          v_surf[i](j)   += (d-d1)/d;
          v_surf[i+1](j) += (d-d2)/d;
          break;
        }
        
      }
    }
  }
  
  // normalize surfaces (like a probability distribution)
  for (int i=0; i<lut_translate.n_slices; i++) {
    for (int j=0; j<n_surface; j++) {
      v_surf[i](j) /= n_pointcloud;
    }
  }
  
  // create shpm (Fourier transfroms) of the surfaces
  VecShpm v_shpm(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    v_shpm[i] = surf2shpm(v_surf[i], contin, l_max);
  }
  
  return v_shpm;

}

