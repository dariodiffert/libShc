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

double Shc::calc_rot_U(MatrixRotation &R, MatrixReal &M, int l, int m, int n) {
  
  double res;
  
  // Fallunterscheidung unnoetig !?
  if (abs(m) > 0) {
    res = calc_rot_P(R,M,0,l,m,n);
  } else {
    res = calc_rot_P(R,M,0,l,0,n);
  }
  
  return res;
  
}
double Shc::calc_rot_V(MatrixRotation &R, MatrixReal &M, int l, int m, int n) {
  
  double res;
  double d1=0; double d2=0;
  
  if (m == -1)
    d1 = 1;
  
  if (m == 1)
    d2 = 1;
  
  
  if (m < 0) {
    res  = calc_rot_P(R,M, 1,l, m+1,n) * (1-d1);
    res += calc_rot_P(R,M,-1,l,-m-1,n) * sqrt(1+d1); // due to symmetries I would have expected a plus sign in the sqrt
  } else if (m > 0) {
    res  = calc_rot_P(R,M, 1,l, m-1,n) * sqrt(1+d2);
    res -= calc_rot_P(R,M,-1,l,-m+1,n) * (1-d2);
  } else {
    res  = calc_rot_P(R,M, 1,l, 1,n);
    res += calc_rot_P(R,M,-1,l,-1,n);
  }
  
  return res;
  
}
double Shc::calc_rot_W(MatrixRotation &R, MatrixReal &M, int l, int m, int n) {
  
  double res=0;
  
  if (m < 0) {
    res  = calc_rot_P(R,M, 1,l, m-1,n);
    res -= calc_rot_P(R,M,-1,l,-m+1,n);
  } else if (m > 0) {
    res  = calc_rot_P(R,M, 1,l, m+1,n);
    res += calc_rot_P(R,M,-1,l,-m-1,n);
  } else {
    print("[ERROR] calc_W: this should not be reachable!\n");
  }
  
  return res;
  
}
double Shc::calc_rot_P(MatrixRotation &R, MatrixReal &M, int i, int l, int a, int b) {
  
  double res = 0;
  
  int x = (l-1)*(l-1)+(l-1);
  
  if (abs(b) < l) {
    res = R(i+1,0+1) * M(a+x,b+x);
  } else if (b == l) {
    res  = R(i+1, 1+1) * M(a+x, l-1+x);
    res -= R(i+1,-1+1) * M(a+x,-l+1+x);
  } else if (b == -l) {
    res  = R(i+1, 1+1) * M(a+x,-l+1+x);
    res += R(i+1,-1+1) * M(a+x, l-1+x);
  } else {
    print("[ERROR] calc_P: this should not be reachable!\n");
  }
  
  return res;
  
}
double Shc::calc_rot_u(int l, int m, int n) { 
  
  double num; double den; double res;
  
  if (abs(n) < l) {
    num = (l+m) * (l-m);
    den = (l+n) * (l-n);
    res = sqrt(num/den);
  } else {
    num = (l+m) * (l-m);
    den = (2*l) * (2*l-1);
    res = sqrt(num/den);
  }
  
  return res;
  
}
double Shc::calc_rot_v(int l, int m, int n) {
  
  double num; double den; double res; double d;
  
  if (m==0) {
    d = 1;
  } else {
    d = 0;
  }
  
  if (abs(n) < l) {
    num = (1+d) * (l+abs(m)-1) * (l+abs(m));
    den = (l+n) * (l-n);
    res = 1.0/2 * (1-2*d) * sqrt(num/den);
  } else {
    num = (1+d) * (l+abs(m)-1) * (l+abs(m));
    den = (2*l) * (2*l-1);
    res = 1.0/2 * (1-2*d) * sqrt(num/den);
  }
  
  return res;
  
}
double Shc::calc_rot_w(int l, int m, int n) {
  
  double num; double den; double res; double d;
  
  if (m==0) {
    d = 1;
  } else {
    d = 0;
  }
  
  if (abs(n) < l) {
    num = (l-abs(m)-1) * (l-abs(m));
    den = (l+n) * (l-n);
    res = -1.0/2.0 * (1-d) * sqrt(num/den);
  } else {
    num = (l-abs(m)-1) * (l-abs(m));
    den = (2*l) * (2*l-1);
    res = -1.0/2.0 * (1-d) * sqrt(num/den);
  }
  
  return res;
  
}
double Shc::calc_rotation_entry(MatrixRotation &R, MatrixReal &M, int l, int m, int n) {
  
  double u = calc_rot_u(l,m,n);
  double v = calc_rot_v(l,m,n);
  double w = calc_rot_w(l,m,n);
   
  double res = 0; 
  
  if (abs(u) > 1e-12) {
    res += u * calc_rot_U(R,M,l,m,n);
  }
    
  if (abs(v) > 1e-12) {
    res += v * calc_rot_V(R,M,l,m,n);
  }
  
  if (abs(w) > 1e-12) {
    res += w * calc_rot_W(R,M,l,m,n);
  }
  
  return res;
  
}
double Shc::calc_rotation_entry_Z(double a,int mi, int mj) {
  
  if (mi * mj >= 0) {
    return cos(abs(mi)*a);
  } else { 
    return sin(mj*a);
  }
  
}

vector<Shc::s_translate_prepare> Shc::calc_translation_prepare(s_coor3 trans, float slice, e_trans_type trans_type) {
   
  s_sphericalCoor sphericalCoor;
  s_coor3 coor;
  float mu, x, y, z, l, u, v, t, s;
  vector<s_translate_prepare> prepare(lut_surf.n_angles);
  
  // for each angle ...
  for (int a=0; a<lut_surf.n_angles; a++) {
    
    // set entries to 'unused'
    prepare[a].mu[0] = 0;
    prepare[a].mu[1] = 0;
    prepare[a].w[0] = 0;
    prepare[a].w[1] = 0;
    
    x = lut_surf.unit[a].coor.x;
    y = lut_surf.unit[a].coor.y;
    z = lut_surf.unit[a].coor.z;
    u = trans.x;
    v = trans.y;
    t = trans.z;
    l = sqrt(u*u+v*v+t*t);
    
    s = u*u*(x*x-1)+v*v*(y*y-1)+t*t*(z*z-1)+2*(u*x*v*y+u*x*t*z+v*y*t*z)+slice*slice;
    
    // ... if there is at least one intersection of the viewing ray with the slice ...
    if (s >= 0) {
      
      // ... calculate for both intersections ...
      for (int i=0; i<2; i++) {
      
        // ... the distances to those.
        mu = -(u*x+v*y+t*z) + pow(-1,i+1) * sqrt(s);
        prepare[a].mu[i] = mu;
        
        // if the new distance is positive (we dont look 'behind') consider this intersection point as valid
        if (mu > 0) {
          
          coor.x = u + mu * x;
          coor.y = v + mu * y;
          coor.z = t + mu * z;

          // calculate the weight of the intersection depending on the distance.
          switch (trans_type) {
            case VISUAL: // only the nearest intersection is used
              if (i==1 && l <= slice) {
                prepare[a].w[i] = 1;
              }
              break;
            case DENSITY: // both intersections are summed up. Weight them by their distance as if they represent a density.
              prepare[a].w[i] = (mu/slice)*(mu/slice);
              break;
          }

          // calculate the surface at the new position (l,0,0) with the given viewing angle
          sphericalCoor = coor32sphericalCoor(coor);
          // apply weighting
          prepare[a].indices[i] = get_nearest_indices(sphericalCoor);
          prepare[a].weights[i] = get_nearest_weights(sphericalCoor, prepare[a].indices[i]);
          
          // if the translation is on the sphere border avoid double summation
          if (s == 0) {
            break;
          }

        }
      }  
      
    }
  }
 
  return prepare;
  
}
vector<Shc::s_translate_prepare> Shc::calc_translation_prepare(float distance, float slice, e_trans_type trans_type) {
   
  s_sphericalCoor sphericalCoor;
  s_coor3 coor;
  float mu, x, y, z, l, s;
  vector<s_translate_prepare> prepare(lut_surf.n_angles);
  
  // for each angle ...
  for (int a=0; a<lut_surf.n_angles; a++) {
    
    // set entries to 'unused'
    prepare[a].mu[0] = 0;
    prepare[a].mu[1] = 0;
    prepare[a].w[0] = 0;
    prepare[a].w[1] = 0;
    
    x = lut_surf.unit[a].coor.x;
    y = lut_surf.unit[a].coor.y;
    z = lut_surf.unit[a].coor.z;
    l = distance;
    s = l*l*(z*z-1)+slice*slice;
    
    // ... if there is at least one intersection of the viewing ray with the slice ...
    if (s >= 0) {
      
      // ... calculate for both intersections ...
      for (int i=0; i<2; i++) {
      
        // ... the distances to those.
        mu = -l*z + pow(-1,i+1) * sqrt(s);
        prepare[a].mu[i] = mu;
        
        // if the new distance is positive (we dont look 'behind') consider this intersection point as valid
        if (mu > 0) {
          
          coor.x = 0 + mu * x;
          coor.y = 0 + mu * y;
          coor.z = l + mu * z;

          // calculate the weight of the intersection depending on the distance.
          switch (trans_type) {
            case VISUAL: // only the nearest intersection is used
              if (i==1 && distance <= slice) {
                prepare[a].w[i] = 1;
              }
              break;
            case DENSITY: // both intersections are summed up. Weight them by their distance as if they represent a density.
              prepare[a].w[i] = (mu/slice)*(mu/slice);
              break;
          }

          // calculate the surface at the new position (l,0,0) with the given viewing angle
          sphericalCoor = coor32sphericalCoor(coor);
          // apply weighting
          prepare[a].indices[i] = get_nearest_indices(sphericalCoor);
          prepare[a].weights[i] = get_nearest_weights(sphericalCoor, prepare[a].indices[i]);
          
          // if the translation is on the sphere border avoid double summation
          if (s == 0) {
            break;
          }

        }
      }  
      
    }
  }
 
  return prepare;
  
}
VecShpm Shc::calc_translation_column(vector<s_translate_prepare>& prepare, int l_max, int index, bool z_only, e_contin contin) {
  
  // create surface of current sh basis function
  Surf surf = coef2surf(index);
  
  // create a surface for each slice, where the translated points are projected on
  VecSurf surf_translate(lut_translate.n_slices);
  for (int b=0; b<lut_translate.n_slices; b++) {
    surf_translate[b].resize(lut_surf.n_angles);
    surf_translate[b].setZero();
  }
  
  float mu, w, val;
  vector<int> indices;
  vector<float> weights;
  
  // for each angle ...
  for (int a=0; a<lut_surf.n_angles; a++) {
    
    // ... calculate for both intersections ...
    for (int i=0; i<2; i++) {
    
      mu = prepare[a].mu[i];
      
      // if the distance is positive (we dont look 'behind') consider this intersection point as valid
      if (mu > 0) {
        
        w       = prepare[a].w[i];
        indices = prepare[a].indices[i];
        weights = prepare[a].weights[i];
        
        // apply weighting
        val = w * get_nearest_surf(surf, indices, weights);
        
        // distribute this new surface information on the surface slices which represent the transformation
        if (mu <= lut_translate.slices(0)) {
          surf_translate[0](a) += val;   
        } else if (mu >= lut_translate.slices(lut_translate.n_slices-1)) {
          surf_translate[lut_translate.n_slices-1](a) += val;   
        } else {

          for (int b=0; b<lut_translate.n_slices-1; b++) {
            
            if (mu >= lut_translate.slices(b) && mu <= lut_translate.slices(b+1)) {
              float d1 = mu - lut_translate.slices(b);
              float d2 = lut_translate.slices(b+1) - mu;
              float d = d1+d2;
              surf_translate[b](a)   += val * (d-d1)/d;
              surf_translate[b+1](a) += val * (d-d2)/d;
              break;
            }
  
          }
        }

      } 

    }
  }
  
  // fourier transform the surfaces back into the spherical harmonics
  VecShpm result(lut_translate.n_slices);
  
  // create transformation matrix
  for (int b=0; b<lut_translate.n_slices; b++) {
    
    result[b] = intern_surf2shpm(surf_translate[b], contin, l_max, false, FULL);
    
    // if set by user: set small values to zero
    int li, mi, lj, mj;
    index2sh(index, li, mi);
    int size = result[b].coef.size();
    for (int i=0; i<size; i++) {
      index2sh(i, lj, mj);
      if (z_only == true && calc_translate_isEntryZero(mi,mj) == true) {
        result[b].coef(i) = 0;
      }      
      if (abs(result[b].coef(i)) < tolerances.translation) {
        result[b].coef(i) = 0;
      }
    }
        
  }

  return result;
  
}
VecMatrixReal Shc::calc_translation(float distance, float slice, e_trans_type trans_type, int l_max) {
  
  VecMatrixReal M(lut_translate.n_slices);
  int sh_max = l2sh(n_bands, FULL);
  for (int i=0; i<lut_translate.n_slices; i++) {
    M[i].resize(sh_max, sh_max);
    M[i].setZero();
  }
  
  VecShpm v_shpm;
  print_progress_start();
  vector<s_translate_prepare> prepare = calc_translation_prepare(distance, slice, trans_type);
  for (int l=0; l<l_max; l++) {
    for (int m=-l; m<=l; m++) {
      int index = sh2index(l, m);
      
      v_shpm = calc_translation_column(prepare, l_max, index, true, FULL);
      for (int i=0; i<lut_translate.n_slices; i++) {
        M[i].col(index) = v_shpm[i].coef;
      }
      print_progress_update(index * 1.0 / sh_max);
    }
  }
  print_progress_stop();
  
  return M;
    
}
MatrixReal Shc::calc_translation_simple(s_coor3 trans, e_trans_type trans_type, e_contin contin, int l_max) {
  
  if (l_max < 0) {
    l_max = n_bands;
  }
  
  int sh_max = l2sh(n_bands, FULL);
  MatrixReal M(sh_max, sh_max);
  M.setZero();
  
  VecShpm v_shpm;
  vector<s_translate_prepare> prepare = calc_translation_prepare(trans, 1, trans_type);
  for (int l=0; l<l_max; l++) {
    for (int m=-l; m<=l; m++) {
      int index = sh2index(l, m);
      
      v_shpm = calc_translation_column(prepare, l_max, index, false, contin);
      M.col(index) = v_shpm[0].coef;
    }
  }
  
  return M;
    
}

MatrixReal Shc::reduce_matrix(MatrixReal& M, e_contin contin) {
  
  if (trifold(contin) == FULL) {
    return M;
  }
  
  int countM; int countN;
  int sh_full = l2sh(n_bands, FULL);
  int sh_cont = l2sh(n_bands, contin);
  
  MatrixReal NC(sh_full, sh_cont);
  NC.setZero();
  
  // Remove all unnecessary columns
  countM = 1;
  countN = 1;

  NC(0,0) = M(0,0);
  for (int l=1; l<n_bands; l++) {;
    
    int n = 2*l+1;
    
    if ((contin == HEMI_RM && l%2 == 0) || (contin == HEMI_RMN && l%2 == 1)) {
      NC.block(0,countN,sh_full,n) = M.block(0,countM,sh_full,n);
      countN += n;
    } 
    countM += n;
  }
 
  return NC;
  
}

int Shc::intern_add_transform(s_transform& transform, MatrixReal& matrix_transform, e_matrix_type matrix_type, float tolerance, int l_max) {
  
  if (l_max < 0) {
    l_max = n_bands;
  }
  
  if (matrix_transform.rows() != n_bands*n_bands || matrix_transform.cols() != n_bands*n_bands) {
    stringstream ss;
    ss << "The passed transformation matrix must have dimension n_bands^2 x n_bands^2 (" << n_bands*n_bands << "x" << n_bands*n_bands << ").";
    print_warning("intern_add_transform", ss.str());
    return -1;
  }
  
  MatrixReal M = matrix_transform;
  
  if (tolerance > 0) {
    int n = M.cols()*M.rows();
    for (int i=0; i<n; i++) {
      if (abs(M(i)) < tolerance) {
        M(i) = 0;
      }
    }
  }
  
  MatrixReal M_RM  = reduce_matrix(M, HEMI_RM);  
  MatrixReal M_RMN = reduce_matrix(M, HEMI_RMN);
  
  s_transform_unit transform_unit;
  transform_unit.l_max       = l_max;
  transform_unit.matrix_type = matrix_type;
  
  transform_unit.M = M;
  
  if (matrix_type == SPARSE) {
    transform_unit.MS     = M.sparseView();
    transform_unit.MS_RM  = M_RM.sparseView();
    transform_unit.MS_RMN = M_RMN.sparseView();
  } else {
    transform_unit.M_RM   = M_RM;
    transform_unit.M_RMN  = M_RMN;
  }
  
  transform.unit.push_back(transform_unit);
  
  return transform.unit.size() - 1;
  
}
void Shc::intern_clear_transform(s_transform& transform) {
  
  transform.unit.clear();
  
}

bool Shc::intern_save_transform(s_transform& transform, string file) {
  
  ofstream fileStream(file, ios::out | ios::binary | ios::trunc);
  
  if (fileStream.is_open()) {
    
    int n = transform.unit.size();
    fileStream.write((char*) (&n), sizeof(int)); 
    for (int i=0; i<n; i++) {
      
      fileStream.write((char*) (&transform.unit[i].l_max),       sizeof(int)); 
      fileStream.write((char*) (&transform.unit[i].matrix_type), sizeof(e_matrix_type)); 
      
      int rows = transform.unit[i].M.rows();
      int cols = transform.unit[i].M.cols();
      fileStream.write((char*) (&rows), sizeof(int));
      fileStream.write((char*) (&cols), sizeof(int));
      fileStream.write((char*) transform.unit[i].M.data(), rows*cols*sizeof(float));

    }
    
    fileStream.close();
    return true;
  }

  return false;
  
}
bool Shc::intern_load_transform(s_transform& transform, string file) {
  
  ifstream fileStream(file, ios::in | ios::binary);
  
  if (fileStream.is_open()) {
    
    int n = 0;
    fileStream.read((char*) (&n), sizeof(int));
    intern_clear_transform(transform);
    
    for (int i=0; i<n; i++) {
      
      s_transform_unit transform_unit;
      
      int l_max;
      e_matrix_type matrix_type;
      fileStream.read((char*) (&l_max),       sizeof(int));
      fileStream.read((char*) (&matrix_type), sizeof(e_matrix_type));
      
      int rows; int cols;
      fileStream.read((char*) (&rows), sizeof(int));
      fileStream.read((char*) (&cols), sizeof(int));
      
      MatrixReal matrix(rows, cols);
      fileStream.read((char*) matrix.data(), rows*cols*sizeof(float));
            
      intern_add_transform(transform, matrix, matrix_type, 0, l_max);
      
    }
    
    fileStream.close();
    return true;
  }
  
  return false;
  
}

Shpm Shc::intern_product(Shpm& shpm1, Shpm& shpm2) {
  
  CoefBandwise bandwise1 = shpm2bandwise(shpm1);
  CoefBandwise bandwise2 = shpm2bandwise(shpm2);
  
  int l1_max = min(n_bands_CG, shpm1.l_max);
  int l2_max = min(n_bands_CG, shpm2.l_max);
  int l_max  = min(n_bands_CG, l1_max+l2_max-1);
  int sh_max = l2sh(l_max, FULL);
 
  Shpm result;
  result.coef.resize(sh_max);
  result.coef.setZero();
  result.l_max  = l_max;
  
  if (shpm1.contin == shpm2.contin) {
    result.contin = shpm1.contin;
  } else {
    result.contin = FULL;
  }
  
  MatrixReal current;
  int i1; int i2;
  
  for (int l1=0; l1<l1_max; l1++) {
    for (int l2=0; l2<l2_max; l2++) {
      
      if (l1 != 0) {
        if (shpm1.contin == HEMI_RM  && odd(l1) )
          continue;
        
        if (shpm1.contin == HEMI_RMN && even(l1) )
          continue;
      }
      
      if (l2 != 0) {
        if (shpm2.contin == HEMI_RM  && odd(l2) )
          continue;
        
        if (shpm2.contin == HEMI_RMN && even(l2) )
          continue;
      }
      
      VectorReal F = kroneckerProduct(bandwise1[l1],bandwise2[l2]);
      
      i1 = (l2-l1)*(l2-l1);
      i2 = (l2+l1+1)*(l2+l1+1)-1;
      
      if (i2 >= sh_max-1)
        i2 = sh_max-1;
      
      int n = i2-i1+1;    
      int m = lut_cg_realCoupling[l1][l2].cols(); // (2*l1+1)*(2*l2+1)

      if (m > n) {
        current = lut_cg_realCoupling[l1][l2].block(0,0,n,n).transpose() * F.segment(0,n) * rescale_factor;
      } else {
        current = lut_cg_realCoupling[l1][l2].transpose() * F * rescale_factor;
      }
      
      result.coef.segment(i1,n) += current.block(0,0,1,n);
      
    }
  }
   
  return result;
  
}
VecShpm Shc::intern_product(VecShpm& v_shpm1, VecShpm& v_shpm2) {

  int n = v_shpm1.size();
  VecShpm v_res(n);
  
  for (int i=0; i<n; i++) {
    v_res[i] = intern_product(v_shpm1[i], v_shpm2[i]);
  }
    
  return v_res;
  
}
VecShpm Shc::intern_rotate(VecShpm& v_shpm, Xyz xyz) {

  s_compass_prepare_shpm scp = rotate_prepare_shpm(v_shpm);
  VecShpm v_result = v_shpm;
  s_rotation_indices indices;

  // find rotation parameterset with highest number of bands for fixed rotations
  int l_max = 0;
  int l_max_index = 0;
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (v_rotation[index_rot_par].l_max > l_max) {
      l_max = v_rotation[index_rot_par].l_max;
      l_max_index = index_rot_par;
    }
  }
  if (l_max == 0) {
    print_warning("intern_rotate", "No rotations (with band L>0) defined");
    return v_result;
  }
  
  // find rotation parameterset with finest resolution
  float resolution_min = -1;
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (resolution_min == -1 || (v_rotation[index_rot_par].resolution != -1 && v_rotation[index_rot_par].resolution < resolution_min)) {
      resolution_min = v_rotation[index_rot_par].resolution;
    }
  }
  
  // rotate vshpm
  float diffz = angular_difference(xyz.z, 0);
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (v_rotation[index_rot_par].l_max > 0) {
      while (diffz > v_rotation[index_rot_par].resolution * 0.51) {
        xyz2rotationIndices(xyz, index_rot_par, indices);
        int i = v_rotation[index_rot_par].fast_indices_z[indices.i];
        if (v_rotation[index_rot_par].psi[i] != 0) {
          rotate_calc_Z(v_result, i, scp, index_rot_par);
          xyz.z -= v_rotation[index_rot_par].psi[i];
          diffz = angular_difference(xyz.z, 0);
        }
        if (v_rotation[index_rot_par].resolution == -1) {break;}
      }
    }
  }

  rotate_calc_ZY(v_result, scp, l_max_index);

  float diffy = angular_difference(xyz.y, 0);
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (v_rotation[index_rot_par].l_max > 0) {
      while (diffy > v_rotation[index_rot_par].resolution * 0.51) {
        xyz2rotationIndices(xyz, index_rot_par, indices);
        int j = v_rotation[index_rot_par].fast_indices_y[indices.j];
        if (v_rotation[index_rot_par].psi[j] != 0) {
          rotate_calc_Z(v_result, j, scp, index_rot_par);
          xyz.y -= v_rotation[index_rot_par].psi[j];
          diffy = angular_difference(xyz.y, 0);
        }
        if (v_rotation[index_rot_par].resolution == -1) {break;}
      }
    }
  }
  
  rotate_calc_YX(v_result, scp, l_max_index);

  float diffx = angular_difference(xyz.x, 0);
  for (int index_rot_par=0; index_rot_par<n_rot_par; index_rot_par++) {
    if (v_rotation[index_rot_par].l_max > 0) {
      while (diffx > v_rotation[index_rot_par].resolution * 0.51) {
        xyz2rotationIndices(xyz, index_rot_par, indices);
        int k = v_rotation[index_rot_par].fast_indices_x[indices.k];
        if (v_rotation[index_rot_par].psi[k] != 0) {
          rotate_calc_Z(v_result, k, scp, index_rot_par);
          xyz.x -= v_rotation[index_rot_par].psi[k];
          diffx = angular_difference(xyz.x, 0);
        }
        if (v_rotation[index_rot_par].resolution == -1) {break;}
      }
    }
  }
  
  rotate_calc_XZ(v_result, scp, l_max_index);

  if (resolution_min != -1) {
    if (diffz > resolution_min || diffy > resolution_min || diffx > resolution_min) {
      stringstream ss;
      ss << "The initialized rotation parameters are not sufficient to perform the desired rotation. ";
      ss << "Did you forget to initialize rotations at all; note that for some functions as translations, rotations need to be initialized! ";
      ss << "(remainer in deg: ";
      ss << "X: " << diffx*180.0f/M_PI << " Y: " <<  diffy*180.0f/M_PI << " Z: " << diffz*180.0f/M_PI << ")";
      print_warning("intern_rotate", ss.str());
    }
  }
  
  return v_result;
  
}
VecShpm Shc::intern_translate(VecShpm& v_shpm, float dist) {
  
  if (lut_translate.n_distances <= 0) {
    print_warning("intern_translate", "No translations defined, use init_translations()!");
    return v_shpm;
  }
  
  VecShpm result;  
  if (lut_translate.n_distances == 1 || bool_linearize_translations == false) {
    result = intern_translate_single(v_shpm, dist);
  } else {
    result = intern_translate_dual(v_shpm, dist);
  }
  
  return result;
  
}
VecShpm Shc::intern_translate_single(VecShpm& v_shpm, float dist) {
  
  float curDist;
  float minDist = 1e12;
  int minIndex = 0;
  for (int i=0; i<lut_translate.n_distances; i++) {
    curDist = abs(lut_translate.distances(i) - dist);
    if (curDist < minDist) {
      minDist = curDist;
      minIndex = i;
    }
  }
  
  vector<int> l_max(lut_translate.n_slices);
  vector<int> sh_max_in(lut_translate.n_slices);
  vector<int> sh_max_res(lut_translate.n_slices);
  
  for (int i=0; i<lut_translate.n_slices; i++) {
    l_max[i]  = min(lut_translate.l_max, v_shpm[i].l_max);
    sh_max_in[i]  = l2sh(l_max[i], v_shpm[i].contin);
    sh_max_res[i] = l2sh(l_max[i], FULL);
  }
  
  VecCoef resultTranslate(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    resultTranslate[i].resize(sh_max_res[i]);
    resultTranslate[i].setZero();
  }
  
  vector<vector<vector<MatrixRealSparse> >*> sparse(lut_translate.n_slices);
  
  for (int j=0; j<lut_translate.n_slices; j++) {
    switch (v_shpm[j].contin) {
      case HEMI_RM:
        sparse[j] = &lut_translate.unit[minIndex].matrixSparse_RM;
        break;
      case HEMI_RMN:
        sparse[j] = &lut_translate.unit[minIndex].matrixSparse_RMN;
        break;
      default:
        sparse[j] = &lut_translate.unit[minIndex].matrixSparse;
        break;
    }
  }
  
  vector<bool> block(lut_translate.n_slices);
  for (int j=0; j<lut_translate.n_slices; j++) {
    if (l2sh(lut_translate.l_max, v_shpm[j].contin) != sh_max_in[j]) {
      block[j] = true;
    } else {
      block[j] = false;
    }
  }
  
  for (int j=0; j<lut_translate.n_slices; j++) {
    for (int k=0; k<lut_translate.n_slices; k++) {
      if (block[j] == true) {
        resultTranslate[k] += ((*(sparse[j]))[j][k]).block(0,0,sh_max_res[j],sh_max_in[j]) * v_shpm[j].coef.segment(0,sh_max_in[j]);
      } else {
        resultTranslate[k] += (*(sparse[j]))[j][k] * v_shpm[j].coef;
      }
    }
  }
  
  VecShpm result(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    result[i].coef   = resultTranslate[i];
    result[i].l_max  = l_max[i];
    result[i].contin = FULL;
  }
  
  return result;
  
}
VecShpm Shc::intern_translate_dual(VecShpm& v_shpm, float dist) {
  
  float curDist;
  float minDist1 = 1e12;
  int minIndex1 = 0;
  for (int i=0; i<lut_translate.n_distances; i++) {
    curDist = abs(lut_translate.distances(i) - dist);
    if (curDist < minDist1) {
      minDist1 = curDist;
      minIndex1 = i;
    }
  }
  
  float minDist2 = 1e12;
  int minIndex2 = 0;
  if (minIndex1 > 0 && minIndex1 < lut_translate.n_distances-1) {
    float t1, t2;
    t1 = abs(lut_translate.distances(minIndex1-1) - dist);
    t2 = abs(lut_translate.distances(minIndex1+1) - dist);
    
    if (t1 < t2) {
      minIndex2 = minIndex1-1;
      minDist2 = t1;
    } else {
      minIndex2 = minIndex1+1;
      minDist2 = t2;
    }
    
  } else if (minIndex1 == 0) {
    minIndex2 = 1;
    minDist2 = abs(lut_translate.distances(minIndex2) - dist);
  } else if (minIndex1 == lut_translate.n_distances-1) {
    minIndex2 = lut_translate.n_distances-2;
    minDist2 = abs(lut_translate.distances(minIndex2) - dist);
  }
  
  float w1 = 1 - minDist1 / (minDist1 + minDist2);
  float w2 = 1 - minDist2 / (minDist1 + minDist2);
  
  vector<int> l_max(lut_translate.n_slices);
  vector<int> sh_max_in(lut_translate.n_slices);
  vector<int> sh_max_res(lut_translate.n_slices);
  
  for (int i=0; i<lut_translate.n_slices; i++) {
    l_max[i]      = min(lut_translate.l_max, v_shpm[i].l_max);
    sh_max_in[i]  = l2sh(l_max[i], v_shpm[i].contin);
    sh_max_res[i] = l2sh(l_max[i], FULL);
  }
  
  VecCoef resultTranslate(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    resultTranslate[i].resize(sh_max_res[i]);
    resultTranslate[i].setZero();
  }
  
  vector<vector<vector<MatrixRealSparse> >*> sparse1(lut_translate.n_slices);
  vector<vector<vector<MatrixRealSparse> >*> sparse2(lut_translate.n_slices);
  
  for (int j=0; j<lut_translate.n_slices; j++) {
    switch (v_shpm[j].contin) {
      case HEMI_RM:
        sparse1[j] = &lut_translate.unit[minIndex1].matrixSparse_RM;
        break;
      case HEMI_RMN:
        sparse1[j] = &lut_translate.unit[minIndex1].matrixSparse_RMN;
        break;
      default:
        sparse1[j] = &lut_translate.unit[minIndex1].matrixSparse;
        break;
    }
  }
  for (int j=0; j<lut_translate.n_slices; j++) {
    switch (v_shpm[j].contin) {
      case HEMI_RM:
        sparse2[j] = &lut_translate.unit[minIndex2].matrixSparse_RM;
        break;
      case HEMI_RMN:
        sparse2[j] = &lut_translate.unit[minIndex2].matrixSparse_RMN;
        break;
      default:
        sparse2[j] = &lut_translate.unit[minIndex2].matrixSparse;
        break;
    }
  }
  
  vector<bool> block(lut_translate.n_slices);
  for (int j=0; j<lut_translate.n_slices; j++) {
    if (l2sh(lut_translate.l_max, v_shpm[j].contin) != sh_max_in[j]) {
      block[j] = true;
    } else {
      block[j] = false;
    }
  }
  
  for (int j=0; j<lut_translate.n_slices; j++) {
    for (int k=0; k<lut_translate.n_slices; k++) {
      if (block[j] == true) {
        resultTranslate[k] += (w1 * ((*(sparse1[j]))[j][k]).block(0,0,sh_max_res[j],sh_max_in[j]) + w2 * ((*(sparse2[j]))[j][k]).block(0,0,sh_max_res[j],sh_max_in[j]) ) * v_shpm[j].coef.segment(0,sh_max_in[j]);
      } else {
        resultTranslate[k] += (w1 * (*(sparse1[j]))[j][k] + w2 * (*(sparse2[j]))[j][k]) * v_shpm[j].coef;
      }
    }
  }
  
  VecShpm result(lut_translate.n_slices);
  for (int i=0; i<lut_translate.n_slices; i++) {
    result[i].coef   = resultTranslate[i];
    result[i].l_max  = l_max[i];
    result[i].contin = FULL;
  }
  
  return result;
  
}
VecShpm Shc::intern_warp(VecShpm& v_shpm, Coor3d translation) {

  VecShpm v_result = v_shpm;
  s_coor3 trans_dir(translation);
  
  s_sphericalCoor sc = coor32sphericalCoor(trans_dir);
  
  Xyz xyzA(0, -sc.theta, -sc.phi);
  Xyz xyzB = transpose(xyzA);

  float dist = translation.norm();
  if (dist == 0) {return v_shpm;}

  v_result = intern_rotate(v_result, xyzA);
  v_result = intern_translate(v_result, dist);
  v_result = intern_rotate(v_result, xyzB);
  
  return v_result;
  
}
Shpm Shc::intern_transform(Shpm& shpm, s_transform_unit& transform) {
  
  int l_max  = min3(shpm.l_max, transform.l_max, n_bands);
  int sh_max_in = l2sh(l_max, shpm.contin);
  int sh_max_res = l2sh(l_max, FULL);
  
  Shpm result;
  result.l_max = l_max;
  result.contin = FULL;
  result.coef.resize(sh_max_res);
  result.coef.setZero();
  
  switch (transform.matrix_type) {
  
    case SPARSE: {
      
      switch (trifold(shpm.contin)) {
        case FULL:
          result.coef.segment(0,sh_max_res) = transform.MS.block(0,0,sh_max_res,sh_max_in)     * shpm.coef.segment(0,sh_max_in);
          break;
        case HEMI_RM:
          result.coef.segment(0,sh_max_res) = transform.MS_RM.block(0,0,sh_max_res,sh_max_in)  * shpm.coef.segment(0,sh_max_in);
          break;
        case HEMI_RMN:
          result.coef.segment(0,sh_max_res) = transform.MS_RMN.block(0,0,sh_max_res,sh_max_in) * shpm.coef.segment(0,sh_max_in);
          break;
        default:
          break;
      }
      break;
    }
      
    case DENSE: {
      
      switch (trifold(shpm.contin)) {
        case FULL:
          result.coef.segment(0,sh_max_res) = transform.M.block(0,0,sh_max_res,sh_max_in)     * shpm.coef.segment(0,sh_max_in);
          break;
        case HEMI_RM:
          result.coef.segment(0,sh_max_res) = transform.M_RM.block(0,0,sh_max_res,sh_max_in)  * shpm.coef.segment(0,sh_max_in);
          break;
        case HEMI_RMN:
          result.coef.segment(0,sh_max_res) = transform.M_RMN.block(0,0,sh_max_res,sh_max_in) * shpm.coef.segment(0,sh_max_in);
          break;
        default:
          break;
      }
      break;
    }  
    
    case BANDWISE: {
      
      int c_full = 0;
      int c_hemi = 0;
      for (int l=0; l<l_max; l++) {
        
        int n = 2*l+1;
        if (contin_skip(shpm.contin, l)) {
          c_full += n;
        } else {
        
          switch (trifold(shpm.contin)) {
            case FULL:
              result.coef.segment(c_full,n) = transform.M.block(c_full,c_full,n,n)      * shpm.coef.segment(c_full,n);
              break;
            case HEMI_RM:
              result.coef.segment(c_full,n) = transform.M_RM.block(c_full,c_hemi,n,n)   * shpm.coef.segment(c_hemi,n);
              break;
            case HEMI_RMN:
              result.coef.segment(c_full,n) = transform.M_RMN.block(c_full,c_hemi,n,n)  * shpm.coef.segment(c_hemi,n);
              break;
            default:
              break;
          }
          c_full += n;
          c_hemi+=n;    
          
        }
        
      }
      break;
    }
    
  }
  return result;
  
}

MatrixReal Shc::intern_warp(MatrixReal& matrix, Coor3d translation, Xyz rotation, e_trans_type trans_type) { 
  
  float w = matrix.cols();
  float h = matrix.rows();
  
  MatrixReal out;
  if (w*h == 0) {
    return out;
  } 
  out.resize(h,w);
  
  s_sphericalCoor sc;
  Coor3d v;
  MatrixRotation R = xyz2r(transpose(rotation));
  
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      
      float t = float(y)/(float)(h) *M_PI;
      float p = (1-float(x)/(float)(w)) * 2.0f*M_PI;
      
      float px = sin(t) * cos(p);
      float py = sin(t) * sin(p);
      float pz = cos(t);
      
      float tx = translation(0);
      float ty = translation(1);
      float tz = translation(2);

      float s0 = (tx*px+ty*py+tz*pz);
      float s = s0*s0 + 1.0f - (tx*tx + ty*ty + tz*tz);
      
      out(y,x) = 0;
      
      if (s >= 0) {
      
        // ... calculate for both intersections ...
        for (int i=0; i<2; i++) {
        
          // ... the distances to those.
          float mu = -s0 + pow(-1,i+1) * sqrt(s);
          float weight = 0;
          
          if (mu > 0) {
              
              v(0) = tx + mu * px;
              v(1) = ty + mu * py;
              v(2) = tz + mu * pz;

              v = R * v;
              
              switch (trans_type) {
                case VISUAL: // only the nearest intersection is used
                  if (i==1 && translation.norm() <= 1) {
                    weight = 1;
                  }
                  break;
                case DENSITY: // both intersections are summed up. Weight them by their distance as if they represent a density.
                  weight = mu*mu;
                  break;
              }

              // calculate the surface at the new position (l,0,0) with the given viewing angle
              sc = coor32sphericalCoor(s_coor3(v));
              
              float j2 = sc.theta*h / M_PI;
              float i2 = w * (1.0f-sc.phi/(2.0f*M_PI));
              
              i2 = fmod(fmod(i2, w) + w, w);
              if (j2 < 0) {
                j2 = 0;
              }
              if (j2 > h-1) {
                j2 = h-1;
              }
              
              out(y,x) += weight * filter_bilinear_mat(matrix,h,w,j2,i2);
              
          }       
          
        }
      }
    }
  }
  
  return out;
  
}

