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

bool Shc::save_lut_translate() {
  
  ofstream file(files.translation, ios::out | ios::binary | ios::trunc);
  
  if (file.is_open()) {
    
    file.write((char*) (&lut_translate.n_distances), sizeof(int));
    file.write((char*) (&lut_translate.n_slices),    sizeof(int));
    file.write((char*) (&lut_translate.l_max),       sizeof(int));
    file.write((char*) (&lut_translate.n_points),    sizeof(int));
    file.write((char*) (&lut_translate.contin),      sizeof(e_contin));
    file.write((char*) (&lut_translate.trans_type),  sizeof(e_trans_type));
    file.write((char*) (&tolerances.translation),    sizeof(float));
    
    file.write((char*) lut_translate.distances.data(), lut_translate.n_distances*sizeof(float));
    file.write((char*) lut_translate.slices.data(), lut_translate.n_slices*sizeof(float));
    
    for (int i=0; i<lut_translate.n_distances; i++) {
      
      file.write((char*) (&lut_translate.distances(i)), sizeof(float));
      
      for (int j=0; j<lut_translate.n_slices; j++) {
        for (int k=0; k<lut_translate.n_slices; k++) {
      
          int rows = lut_translate.unit[i].matrixDense[j][k].rows();
          int cols = lut_translate.unit[i].matrixDense[j][k].cols();
          file.write((char*) (&rows), sizeof(int));
          file.write((char*) (&cols), sizeof(int));
          file.write((char*) lut_translate.unit[i].matrixDense[j][k].data(), rows*cols*sizeof(float));
        
        }
      }
    }
    
    file.close();
    return true;
  }
  
  return false;
  
}
bool Shc::load_lut_translate() {
  
  ifstream file(files.translation, ios::in | ios::binary);
  
  if (file.is_open()) {
    
    file.read((char*) (&lut_translate.n_distances), sizeof(int));
    file.read((char*) (&lut_translate.n_slices),    sizeof(int));
    file.read((char*) (&lut_translate.l_max),       sizeof(int));
    file.read((char*) (&lut_translate.n_points),    sizeof(int));
    file.read((char*) (&lut_translate.contin),      sizeof(e_contin));
    file.read((char*) (&lut_translate.trans_type),  sizeof(e_trans_type));
    file.read((char*) (&tolerances.translation),    sizeof(float));
    
    lut_translate.distances.resize(lut_translate.n_distances);
    lut_translate.unit.resize(lut_translate.n_distances);
    
    file.read((char*) lut_translate.distances.data(), lut_translate.n_distances*sizeof(float));
    file.read((char*) lut_translate.slices.data(),    lut_translate.n_slices*sizeof(float));
    
    for (int i=0; i<lut_translate.n_distances; i++) {
      
      file.read((char*) (&lut_translate.distances(i)), sizeof(float));
      
      for (int j=0; j<lut_translate.n_slices; j++) {
        lut_translate.unit[i].matrixDense.resize(lut_translate.n_slices);
        lut_translate.unit[i].matrixSparse.resize(lut_translate.n_slices);
        lut_translate.unit[i].matrixSparse_RM.resize(lut_translate.n_slices);
        lut_translate.unit[i].matrixSparse_RMN.resize(lut_translate.n_slices);
        for (int k=0; k<lut_translate.n_slices; k++) {
          lut_translate.unit[i].matrixDense[j].resize(lut_translate.n_slices);
          lut_translate.unit[i].matrixSparse[j].resize(lut_translate.n_slices);
          lut_translate.unit[i].matrixSparse_RM[j].resize(lut_translate.n_slices);
          lut_translate.unit[i].matrixSparse_RMN[j].resize(lut_translate.n_slices);
          
          int rows; int cols;
          file.read((char*) (&rows), sizeof(int));
          file.read((char*) (&cols), sizeof(int));
          lut_translate.unit[i].matrixDense[j][k].resize(rows, cols);
          file.read((char*) lut_translate.unit[i].matrixDense[j][k].data(), rows*cols*sizeof(float));
          
          MatrixReal M_reduced_RM  = reduce_matrix(lut_translate.unit[i].matrixDense[j][k], HEMI_RM);
          MatrixReal M_reduced_RMN = reduce_matrix(lut_translate.unit[i].matrixDense[j][k], HEMI_RMN);
          
          lut_translate.unit[i].matrixSparse[j][k]     = lut_translate.unit[i].matrixDense[j][k].sparseView();
          lut_translate.unit[i].matrixSparse_RM[j][k]  = M_reduced_RM.sparseView();
          lut_translate.unit[i].matrixSparse_RMN[j][k] = M_reduced_RMN.sparseView();
          
        }
      }
    }
    
    file.close();
    return true;
  }
  
  return false;
  
}
bool Shc::save_lut_surf() {
  
  ofstream file(files.surface, ios::out | ios::binary | ios::trunc);
  
  if (file.is_open()) {
    
    file.write((char*) (&lut_surf.l_max),                sizeof(int));
    file.write((char*) (&lut_surf.n_angles),             sizeof(int));
    file.write((char*) (&lut_surf.n_sha),                sizeof(int));
    file.write((char*) (&lut_surf.type),                 sizeof(e_surf));
    file.write((char*) lut_surf.sh.data(),               lut_surf.n_sha*sizeof(float));
    file.write((char*) lut_surf.sh_weighted_hemi.data(), lut_surf.n_sha*sizeof(float));
    file.write((char*) lut_surf.sh_weighted_full.data(), lut_surf.n_sha*sizeof(float));
    
    for (int i=0; i<lut_surf.n_angles; i++) {
      file.write((char*) (&lut_surf.unit[i].sphericalCoor), sizeof(s_sphericalCoor));
      file.write((char*) (&lut_surf.unit[i].coor),          sizeof(s_coor3));
      file.write((char*) (&lut_surf.unit[i].weight),        sizeof(float));
    }
    
    file.close();
    return true;
  }
  
  return false;
  
}
bool Shc::load_lut_surf() {
  
  ifstream file(files.surface, ios::in | ios::binary);
  
  if (file.is_open()) {
    
    file.read((char*) (&lut_surf.l_max),              sizeof(int));
    file.read((char*) (&lut_surf.n_angles),           sizeof(int));
    file.read((char*) (&lut_surf.n_sha),              sizeof(int));
    file.read((char*) (&lut_surf.type),               sizeof(e_surf));
    
    lut_surf.sh.resize(lut_surf.n_sha);
    lut_surf.sh_weighted_hemi.resize(lut_surf.n_sha);
    lut_surf.sh_weighted_full.resize(lut_surf.n_sha);
    file.read((char*) lut_surf.sh.data(),               lut_surf.n_sha*sizeof(float));
    file.read((char*) lut_surf.sh_weighted_hemi.data(), lut_surf.n_sha*sizeof(float));
    file.read((char*) lut_surf.sh_weighted_full.data(), lut_surf.n_sha*sizeof(float));
    
    lut_surf.unit.resize(lut_surf.n_angles);
    for (int i=0; i<lut_surf.n_angles; i++) {
      file.read((char*) (&lut_surf.unit[i].sphericalCoor), sizeof(s_sphericalCoor));
      file.read((char*) (&lut_surf.unit[i].coor),          sizeof(s_coor3));
      file.read((char*) (&lut_surf.unit[i].weight),        sizeof(float));
    }
    
    file.close();
    return true;
  }
  
  return false;
  
}

vector<Coor3d> Shc::read_file_vec3(string file_name) {
  
  vector<Coor3d> data;

  ifstream file(file_name);

  if (file.is_open() == false) {
    print_warning("read_file_vec3","Could not open file");
    return data;
  }

  while (file) {
    
    string line;
    
    if (!getline(file, line))
      break;
    
    stringstream ss(line);
    Coor3d v;
    
    getline(ss, line, ',');
    stringstream(line) >> v(0);
    getline(ss, line, ',');
    stringstream(line) >> v(1);
    getline(ss, line, ',');
    stringstream(line) >> v(2);
    
    data.push_back(v);

  }
  
  return data;
  
}
