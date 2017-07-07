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

#include <Shx.h>

using namespace std; 
using namespace Eigen;
using namespace shc;

// These functions will load/save/warp/copy images using the Shc class.
// These do not need to be edited/adapted.
MatrixReal Shx::load_matrix(string file) {
 
  stringstream ss;
  ss << "load_matrix: '" << file << "'\n";
  print(ss.str());
  
  MatrixReal matrix;
  
  Image image;
  if(read_image(image, file) == false) {
    print_warning("load_matrix", "Could not read file.");
    return matrix;
  }

  matrix = image2matrix(image);

  return matrix;
  
}
bool Shx::save_matrix(MatrixReal& matrix, string file) {

  stringstream ss;
  ss << "save_matrix: '" << file << "'\n";
  print(ss.str());
  
  Image image = matrix2image(matrix);
  
  bool result = write_image(image, file);
  
  if (result == false) {
    print_warning("save_matrix", "Could not write matrix to file.");
  }
  
  return result;

}
Shpm Shx::load_shpm(string file, e_contin contin, int l_max) {
 
  stringstream ss;
  ss << "load_shpm: '" << file << "'\n";
  print(ss.str());
  
  MatrixReal matrix = load_matrix(file);
  Shpm shpm         = matrix2shpm(matrix, contin, l_max);

  return shpm;
  
}
bool Shx::save_shpm(Shpm& shpm, string file) {
  
  stringstream ss;
  ss << "save_shpm: '" << file << "'\n";
  print(ss.str());
  
  MatrixReal matrix = shpm2matrix(shpm); 
  
  bool result = save_matrix(matrix, file);

  return result;

}
Shpm Shx::load_shpm_ocamcalib(string file, e_contin contin, int l_max) {
 
  stringstream ss;
  ss << "load_shpm_ocamcalib: '" << file << "'\n";
  print(ss.str());
  
  MatrixReal matrix = load_matrix(file);
  Shpm shpm         = ocamcalib2shpm(matrix, contin, l_max);

  return shpm;
  
}

void Shx::show(Shpm& shpm, string title, bool waitKey) {

  MatrixReal matrix = shpm2matrix(shpm);
  show(matrix, title, waitKey);
  
}
void Shx::show(MatrixReal& matrix, string title, bool waitKey) {

  if (matrix.cols() <= 0 || matrix.rows() <= 0) {
    print_warning("show", "matrix is empty.");
    return;
  }
  
  Image image = matrix2image(matrix);
  show(image, title, waitKey);
  
}
void Shx::show(initializer_list<Shpm> shpm, string title, bool waitKey) {
  
  vector<Shpm> v_shpm;
  for (Shpm current : shpm) {v_shpm.push_back(current);}
  show(v_shpm, title, waitKey);
  
}
void Shx::show(initializer_list<MatrixReal> matrix, string title, bool waitKey) {
  
  vector<MatrixReal> v_matrix;
  for (MatrixReal current : matrix) {v_matrix.push_back(current);}
  show(v_matrix, title, waitKey);
  
}
void Shx::show(vector<Shpm>& shpm, string title, bool waitKey) {
  
  int n = shpm.size();
  
  vector<MatrixReal> v_matrix(n);
  
  for (int i=0; i<n; i++) {
    v_matrix[i] = shpm2matrix(shpm[i]);
  }
  
  show(v_matrix, title, waitKey);
  
}
void Shx::show(vector<MatrixReal>& matrix, string title, bool waitKey) {
  
  int n = matrix.size();
  
  if (n == 0) {
    print_warning("show", "An empty list has been passed");
    return;
  }
  
  int h = matrix[0].rows();
  int w = matrix[0].cols();
  
  for (int i=1; i<n; i++) {
    if (matrix[i].rows() != h || matrix[i].cols() != w) {
      print_warning("show", "All matrices in the list need to be of the same dimension");
      return;
    }
  }
  
  const int margin = 3;
  
  MatrixReal res(h, w*n + margin * (n-1));
  res.setZero();
  
  for (int i=0; i<n; i++) {
    res.block(0,(w+margin)*i,h,w) = matrix[i];
  }
  
  show(res, title, waitKey);
  
}

Image Shx::load_image(string file) {
  
  stringstream ss;
  ss << "load_image: '" << file << "'\n";
  print(ss.str());
  
  Image image;
  if(read_image(image, file) == false) {
    print_warning("load_image", "Could not read file.");
  }
  
  return image;
}
bool Shx::save_image(Image& image, string file) {
  
  stringstream ss;
  ss << "save_image: '" << file << "'\n";
  print(ss.str());
  
  bool result = write_image(image, file);
  
  if (result == false) {
    print_warning("write_image", "Could not write image to file.");
  }
  
  return result;
  
}

