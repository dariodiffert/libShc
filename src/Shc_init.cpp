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

bool Shc::init_lib_translations_load() {
  
  s_translate save_translate = lut_translate;
  float save_tolerance_translation = tolerances.translation;
  
  stringstream ss;
  
  // load translation data from file
  if (files.translation != "") {
    if (load_lut_translate() == true) {
      
      print("Successfully loaded data from file. ");
      
      bool differ_slices    = false;
      if (save_translate.n_slices == lut_translate.n_slices) {
        for (int i=0; i<save_translate.n_slices; i++) {
          if (save_translate.slices(i) != lut_translate.slices(i)) {
            differ_slices = true;
            break;
          }
        }
      } else {
        differ_slices = true;
      }
      
      bool differ_distances = false;
      if (save_translate.n_distances == lut_translate.n_distances) {
        for (int i=0; i<save_translate.n_distances; i++) {
          if (save_translate.distances(i) != lut_translate.distances(i)) {
            differ_distances = true;
            break;
          }
        }
      } else {
        differ_distances = true;
      }
      
      if (save_translate.l_max != lut_translate.l_max || differ_distances || differ_slices || lut_surf.n_angles != lut_translate.n_points || save_tolerance_translation != tolerances.translation || save_translate.trans_type != lut_translate.trans_type) {
        
        ss << "Loaded values (l_max,n_step,n_slices,n_points,type,tolerance) are different from set values, consider recalculating (settings ";
        ss << save_translate.l_max << " " << save_translate.n_distances << " " << save_translate.n_slices << " " <<   lut_surf.n_angles   << " " << save_translate.trans_type   << " " << save_tolerance_translation << " loaded ";
        ss <<  lut_translate.l_max << " " <<  lut_translate.n_distances << " " <<  lut_translate.n_slices << " " << lut_translate.n_points   << " " <<  lut_translate.trans_type   << " " <<     tolerances.translation << "). ";
        
        if (differ_distances) {
          ss << "The set translation distance steps differ from the loaded ones. ";
        }
        if (differ_slices) {
          ss << "The set slices radii differ from the loaded ones. ";
        }
        
        if (files.force_overwrite == false) {
          print_warning("init_lib_translations", ss.str());
          return true;
        } else if (bool_sphere_mode == false) {
          ss << "Using loaded file (force_overwrite == true is ignored) since new translation data can only be calculated using a sphere instead of custom sampling.";
          print_warning("init_lib_translations", ss.str());
          return true;
        } else {
          ss << "Values will be recalculated (force_overwrite == true).";
          print_warning("init_lib_translations", ss.str());
          ss.str("");
        }

      } else {
        return true;
      }

    }
  }
  
  tolerances.translation = save_tolerance_translation;
  lut_translate = save_translate;
  
  return false;
  
}
bool Shc::init_lib_surface_load() {
  
  s_surf save_surface = lut_surf;
  
  int temp_l_max = 0;
  
  stringstream ss;
  
  if (files.surface != "") {
    if (load_lut_surf() == true) {
      
      temp_l_max = lut_surf.l_max;
      
      if (lut_surf.l_max < n_bands) {
        lut_surf.l_max    = n_bands;
        lut_surf.n_sha    = l2sh(n_bands, FULL)*lut_surf.n_angles;
        VectorReal store = lut_surf.sh;
        lut_surf.sh.resize(lut_surf.n_sha);
        lut_surf.sh.setZero();
        lut_surf.sh.segment(0, store.size()) = store;
      }
      
      print("Successfully loaded data from file. ");
      if (n_bands != temp_l_max || save_surface.n_angles != lut_surf.n_angles || save_surface.type != lut_surf.type) {
        
        string type1, type2;
        switch (save_surface.type) {
          case SURF_UNIFORM_SPHERE:
            type1 = "\"uniform_sphere\"";
            break;
          case SURF_UNIFORM_PANORAMA:
            type1 = "\"uniform_panorama\"";
            break;
          default:
            type1 = "\"custom\"";
            break;
        }
        switch (lut_surf.type) {
          case SURF_UNIFORM_SPHERE:
            type2 = "\"uniform_sphere\"";
            break;
          case SURF_UNIFORM_PANORAMA:
            type2 = "\"uniform_panorama\"";
            break;
          default:
            type2 = "\"custom\"";
            break;
        }

        ss << "Loaded values (n_bands,n_points,type) are different from set values, consider recalculating (settings ";
        ss <<            n_bands << " " << save_surface.n_angles << " " << type1 << " loaded ";
        ss <<         temp_l_max << " " << lut_surf.n_angles  << " " << type2  << "). ";
        
        if (files.force_overwrite == false) {
          print_warning("init_lib_surface", ss.str());
          return true;
        } else {
          ss << "Values will be recalculated (force_overwrite == true).";
          print_warning("init_lib_surface", ss.str());
          ss.str("");
        }
        
      } else {
        return true;
      }
    }
  }
  
  lut_surf = save_surface;
  
  return false;
  
}
void Shc::init_lib_parameters() {
  
  n_rot_par              = v_rotation.size();
  
  n_rotate_psi.resize(n_rot_par);
  n_rotate_theta.resize(n_rot_par);
  n_rotate_phi.resize(n_rot_par);
    
  for (int q=0; q<n_rot_par; q++) {
    if (v_rotation[q].psi.size() + v_rotation[q].theta.size() + v_rotation[q].phi.size() == 0) {
      v_rotation[q].n_fast_indices_x    = 0;
      v_rotation[q].n_fast_indices_y    = 0;
      v_rotation[q].n_fast_indices_z    = 0;
      v_rotation[q].n_total_rotations   = 0;
      v_rotation[q].l_max               = 0;
      v_rotation[q].l_max_compass       = 0;
    } 
  }
    
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  
  int sh_full = l2sh(n_bands,FULL);
  
  // ----------------------------------- n_sha -----------------------------------
  lut_surf.n_sha           = sh_full*lut_surf.n_angles;
  lut_surf.l_max           = n_bands;
  
  // ----------------------------------- lut_index2sh -----------------------------------
  lut_index2sh.resize(sh_full);
  
  int count_sh = 0;
  for (int l=0; l<n_bands; l++) {
    for (int m=-l; m<=l; m++) {
      lut_index2sh[count_sh].l = l;
      lut_index2sh[count_sh].m = m; 
      count_sh++;
    }
  }
  
}
void Shc::init_lib_output() {
  
  stringstream ss;
  
  int n = n_bands*n_bands*output_height;
  
  lut_IFFT_output.cfg.init(output_width, true);
  lut_IFFT_output.in.resize(output_width);
  lut_IFFT_output.val.resize(n);

  print_progress_start();
  int count = 0;
  for (int l=0; l<n_bands; l++) {
    for (int m=-l; m<=l; m++) {
      for (int h=0; h<output_height; h++) {
        
        float t = (float)(h) / (float)(output_height-1) * M_PI;
        
        if (m==0) {
          lut_IFFT_output.val[count] =     2.0f * M_PI * calc_sh_K(l, 0) * calc_sh_P(l, 0, cos(t));
        } else if (m < 0) {
          lut_IFFT_output.val[count] = -sqrt(2) * M_PI * calc_sh_K(l,-m) * calc_sh_P(l,-m, cos(t));
        } else if (m > 0) {
          lut_IFFT_output.val[count] =  sqrt(2) * M_PI * calc_sh_K(l, m) * calc_sh_P(l, m, cos(t));
        }
        
        print_progress_update(count * 1.0 / n);
        count++;
        
      }
    }
  }
  print_progress_stop();
  
}
void Shc::init_lib_rotations(int index_rot_par) {
  
  n_rotate_psi[index_rot_par]   = v_rotation[index_rot_par].psi.size();
  n_rotate_theta[index_rot_par] = v_rotation[index_rot_par].theta.size();
  n_rotate_phi[index_rot_par]   = v_rotation[index_rot_par].phi.size();
    
  v_rotation[index_rot_par].rotate_Z.resize(n_rotate_psi[index_rot_par]);
    
  if (v_rotation[index_rot_par].l_max < 0) {
    v_rotation[index_rot_par].l_max = n_bands;
  }
  if (v_rotation[index_rot_par].l_max > n_bands) {
    print_warning("init_lib_rotations", "rot_par set to a higher band than n_bands. rot_par will be adjusted");
  }
  
  if (v_rotation[index_rot_par].l_max_compass > v_rotation[index_rot_par].l_max) {
    print_warning("init_lib_rotations", "The value of l_max_compass has to be smaller or equal to l_max.");
    v_rotation[index_rot_par].l_max_compass = v_rotation[index_rot_par].l_max;
  }
  if (v_rotation[index_rot_par].l_max_compass < 0 ) {
    v_rotation[index_rot_par].l_max_compass = v_rotation[index_rot_par].l_max;
  }
      
  int count = 0;
  int count_total = n_rotate_psi[index_rot_par] + n_rotate_theta[index_rot_par] + n_rotate_phi[index_rot_par];
  
  stringstream ss;
  ss << index_rot_par+1 << "/" << n_rot_par << " ";
  print(ss.str());
  ss.str("");
  
  print_progress_start();
  MatrixReal M, MX, MY, MP, MN;
  MatrixRealSparse MXS, MYS, MS;
  MatrixReal MR;
  
  for (int i=0; i<n_rotate_psi[index_rot_par]; i++) {
    v_rotation[index_rot_par].rotate_Z[i] = create_rotate_Z(v_rotation[index_rot_par].psi[i], v_rotation[index_rot_par].l_max);
    print_progress_update(count++ * 1.0 / count_total);
  }
  
  MN = calc_rotate_axis(-M_PI/2.0, AXIS_X, v_rotation[index_rot_par].l_max);
  MP = MN.transpose();
  v_rotation[index_rot_par].rotate_XN = create_rotate_X(MN, v_rotation[index_rot_par].l_max);
  v_rotation[index_rot_par].rotate_XP = create_rotate_X(MP, v_rotation[index_rot_par].l_max);
  
//   s_rot_X2 rot = create_rotate_X2(MN, v_rotation[index_rot_par].l_max);
//   Coef coef(MN.rows()); for (int i=0; i<MN.rows(); i++) {coef(i) = i;}
//   Coef res2 = compute_rotate_X2(rot, coef, v_rotation[index_rot_par].l_max);

  MN = calc_rotate_axis(-M_PI/2.0, AXIS_Y, v_rotation[index_rot_par].l_max);
  MP = MN.transpose();
  v_rotation[index_rot_par].rotate_YN = create_rotate_Y(MN, v_rotation[index_rot_par].l_max);
  v_rotation[index_rot_par].rotate_YP = create_rotate_Y(MP, v_rotation[index_rot_par].l_max);
  
  if (count_total != 0) {
    MX = calc_rotate_axis(-M_PI/2.0, AXIS_X, v_rotation[index_rot_par].l_max);
    MY = calc_rotate_axis(-M_PI/2.0, AXIS_Y, v_rotation[index_rot_par].l_max);
    MXS = MX.sparseView();
    MYS = MY.sparseView();
    MS = MYS * MXS;
    M = MatrixReal(MS);
  }
  
  v_rotation[index_rot_par].rotate_YNXN = create_rotate_X(M, v_rotation[index_rot_par].l_max);
  
  print_progress_stop();
  
}
void Shc::init_lib_tm() {
  
  lut_tm.resize(2*n_bands);
  
  for (int l=0; l<2*n_bands; l++) {
    lut_tm[l] = createTransformationMatrix(l);
  }
  
}
void Shc::init_lib_cg() {
  
  int count = 0;
  int count_total = 0;
  
  lut_cg.resize(n_bands_CG);
  lut_cg_real.resize(n_bands_CG);
  lut_cg_realCoupling.resize(n_bands_CG);
  
  for (int l=0; l<n_bands_CG; l++) {
    lut_cg[l].resize(n_bands_CG);
    lut_cg_real[l].resize(n_bands_CG);
    count_total += n_bands_CG;
    lut_cg_realCoupling[l].resize(n_bands_CG);
  }
  
  print_progress_start();
  for (int k=0; k<n_bands_CG; k++) {
    for (int l=0; l<n_bands_CG; l++) {
      lut_cg[k][l]          = cg_createMatrix(k,l);
      lut_cg_real[k][l]     = cg_createMatrixReal(k,l);
      lut_cg_realCoupling[k][l] = cg_createMatrixRealCoupling(k,l);
      print_progress_update(count++ * 1.0 / count_total);
    }
  }
  print_progress_stop();
  
}
void Shc::init_lib_translations() {
  
  stringstream ss;
  
  if (lut_translate.l_max < 0) {
    lut_translate.l_max = n_bands;
  }
  
  if (lut_translate.l_max > n_bands) {
    lut_translate.l_max = n_bands;
    print_warning("init_lib_translations", "Number of bands for translations is higher than the maximal number of bands and will be adjusted.");
  }
  
  if (lut_translate.n_distances > 0) {
    
    // get translation settings
    lut_translate.n_points = lut_surf.n_angles; 
    
    // create data if loading failed
    if (init_lib_translations_load() == false) {
    
      print("Calculating... ");
      
      if (bool_sphere_mode == false) {
        lut_translate.n_distances = 0;
        print_warning("init_lib_translations", "Translations can only be calculated using non-custom sample points. Precalculate them on sphere first to do so");
        return;
      }
      
      s_translate_unit unit;
      for (int i=0; i<lut_translate.n_distances; i++) {
        unit.matrixDense.resize(lut_translate.n_slices);
        unit.matrixSparse.resize(lut_translate.n_slices);
        unit.matrixSparse_RM.resize(lut_translate.n_slices);
        unit.matrixSparse_RMN.resize(lut_translate.n_slices);
        for (int j=0; j<lut_translate.n_slices; j++) {
          unit.matrixDense[j].resize(lut_translate.n_slices);
          unit.matrixSparse[j].resize(lut_translate.n_slices);
          unit.matrixSparse_RM[j].resize(lut_translate.n_slices);
          unit.matrixSparse_RMN[j].resize(lut_translate.n_slices);
        }
        this->lut_translate.unit.push_back(unit);
      }
        
      for (int i=0; i<lut_translate.n_distances; i++) {
        for (int j=0; j<lut_translate.n_slices; j++) {
          
          ss << i*lut_translate.n_slices+j+1 << "/" << lut_translate.n_distances*lut_translate.n_slices << " ";
          print(ss.str());
          ss.str("");
          
          VecMatrixReal M = calc_translation(lut_translate.distances(i), lut_translate.slices(j), lut_translate.trans_type, lut_translate.l_max);
          VecMatrixReal M_reduced_RM(lut_translate.n_slices);
          VecMatrixReal M_reduced_RMN(lut_translate.n_slices);
          for (int k=0; k<lut_translate.n_slices; k++) {
            
            M_reduced_RM[k]  = reduce_matrix(M[k], HEMI_RM);
            M_reduced_RMN[k] = reduce_matrix(M[k], HEMI_RMN);
            
            lut_translate.unit[i].matrixDense[j][k]      = M[k];
            lut_translate.unit[i].matrixSparse[j][k]     = M[k].sparseView();
            lut_translate.unit[i].matrixSparse_RM[j][k]  = M_reduced_RM[k].sparseView();
            lut_translate.unit[i].matrixSparse_RMN[j][k] = M_reduced_RMN[k].sparseView();

          }
          
        }
      }
      
      // save translation data to file
      if (files.translation != "" && files.force_overwrite == true) {
        if (save_lut_translate() == true) {
          print("Successfully saved translation data to file. ");
        } else {
          print_warning("init_lib_translations", "Could not save translation data to file");
        }
      }
      
    }
    
  }
  
}
void Shc::init_lib_surface_sh(s_surf& surface) { 
  
  // sh allocation  
  surface.sh.resize(surface.n_sha);

  print_progress_start();
  double y;
  int count_sha = 0;
  for (int l=0; l<surface.l_max; l++) {
    for (int m=-l; m<=l; m++) {
      
      for (int a=0; a<surface.n_angles; a++) {
        y = calc_sh_Y(l, m, surface.unit[a].sphericalCoor.theta, surface.unit[a].sphericalCoor.phi);
        if (std::isnan(y)) {
          surface.sh[count_sha] = 0;
        } else {
          surface.sh[count_sha] = y;
        }
        count_sha++;
        if (count_sha % 100 == 0) {
          print_progress_update(count_sha * 1.0 / surface.n_sha);
        }
      }
      
    }
  }
  
  surface.sh_weighted_full.resize(surface.n_sha);
  surface.sh_weighted_hemi.resize(surface.n_sha);
  
  for (int l=0; l<surface.l_max; l++) {
    for (int m=-l; m<=l; m++) {
      for (int a=0; a<surface.n_angles; a++) {
        
        int i = sha2index(surface,l,m,a);
        
        surface.sh_weighted_full[i] = surface.sh[i] * surface.unit[a].weight;
        
        if (surface.unit[a].coor.z < -1e-4) {
          surface.sh_weighted_hemi[i] = 0;
        } else if (surface.unit[a].coor.z < 1e-4) {
          surface.sh_weighted_hemi[i] = surface.sh[i] * surface.unit[a].weight; // no factor 2.0, since on the equator
        } else {
          surface.sh_weighted_hemi[i] = surface.sh[i] * surface.unit[a].weight * 2.0; // factor 2.0 is due to the doubled area (original+mirrored) on the surface, covered by this coefficient
        }
        
      }
    }
  }
  
  
  print_progress_stop();
   
}
void Shc::init_lib_surface_FFT() {

  lut_FFT.theta_index.resize(lut_surf_quickRef.n_theta);
  
  vector<int> fft_n;
  
  int last_fft_n = -1;
  for (int t=0; t<=(lut_surf_quickRef.n_theta-1)/2; t++) {
    int cur_fft_n = lut_surf_quickRef.psi_count[t];
    if (cur_fft_n != last_fft_n) {
      fft_n.push_back(cur_fft_n);
      last_fft_n = cur_fft_n;
    }
    lut_FFT.theta_index[t]                             = fft_n.size()-1;
    lut_FFT.theta_index[lut_surf_quickRef.n_theta-1-t] = fft_n.size()-1;
  }
  
  lut_FFT.cfg_FFT.resize(fft_n.size());
  lut_FFT.cfg_IFFT.resize(fft_n.size());
  int fft_n_n = fft_n.size();
  for (int i=0; i<fft_n_n; i++) {
    lut_FFT.cfg_FFT[i].init(fft_n[i], false);
    lut_FFT.cfg_IFFT[i].init(fft_n[i], true);    
  } 
  
  lut_FFT.store.resize(lut_surf_quickRef.n_theta);
  for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
    int mm = max(lut_surf_quickRef.psi_count[t], n_bands*n_bands);
    lut_FFT.store[t].resize(mm);
  }
  
  int n = n_bands*n_bands*lut_surf_quickRef.n_theta;
  lut_FFT.val_FFT.resize(n);
  lut_FFT.val_IFFT.resize(n);
  
  float weight_FFT = 1;
  float weight_IFFT = 1;
  float val = 1;
  
  print_progress_start();
  int count = 0;
  for (int l=0; l<n_bands; l++) {
    for (int m=-l; m<=l; m++) {
      for (int t=0; t<lut_surf_quickRef.n_theta; t++) {

        if (lut_surf.type == SURF_UNIFORM_PANORAMA) {
          weight_FFT = M_PI / (float)lut_surf.n_angles * sin(lut_surf_quickRef.theta[t]);
        } else {
          weight_FFT = rescale_factor / (float)lut_surf.n_angles;
        }
        
        if (m==0) {
          weight_IFFT = 2.0f * M_PI;
        } else {
          weight_IFFT = M_PI;
        }
        
        if (m==0) {
          val =            calc_sh_K(l, 0) * calc_sh_P(l, 0, cos(lut_surf_quickRef.theta[t]));
        } else if (m < 0) {
          val = -sqrt(2) * calc_sh_K(l,-m) * calc_sh_P(l,-m, cos(lut_surf_quickRef.theta[t]));
        } else if (m > 0) {
          val =  sqrt(2) * calc_sh_K(l, m) * calc_sh_P(l, m, cos(lut_surf_quickRef.theta[t]));
        }
        
        lut_FFT.val_FFT[count]  = weight_FFT  * val;
        lut_FFT.val_IFFT[count] = weight_IFFT * val;
        
        print_progress_update(count * 1.0f / n);
        count++;
        
      }
    }
  }
  print_progress_stop();
  
}
void Shc::init_lib_surface() {
 
  if (bool_sphere_mode == true) {
    lut_surf_quickRef = create_surface_quickRef(lut_surf);
  }
    
  if (bool_FFT == true) {
    init_lib_surface_FFT();
  } else {
    lut_surf.l_max    = n_bands;
    lut_surf.n_angles = lut_surf.unit.size();
    lut_surf.n_sha    = l2sh(n_bands,FULL)*lut_surf.n_angles;
  
    if (init_lib_surface_load() == false) {
      print("Calculating... ");
      init_lib_surface_sh(lut_surf);
      
      if (files.surface != "") {
        if (save_lut_surf() == true) {
          print("Successfully saved surface data to file. ");
        } else {
          print_warning("init_lib_surface", "Could not save surface data to file");
        }
      }
      
    }
  }
  
  if (bool_init_ocamcalib == true) {
    // create mappings
    ocam_model.mapped.resize(lut_surf.n_angles);   
    for (int i=0; i<lut_surf.n_angles; i++) {
      s_coor2 mapped;
      world2cam(mapped, lut_surf.unit[i].coor);
      ocam_model.mapped[i].x = mapped.y;
      ocam_model.mapped[i].y = mapped.x;
    }
    
    ocam_model.mappingConversionY.resize(output_height, output_width);
    ocam_model.mappingConversionX.resize(output_height, output_width);
    s_coor2 coor2;
    for (float y=0; y<output_height; y++) {
      for (float x=0; x<output_width; x++) {
        
        s_sphericalCoor sc;
        sc.theta = y/output_height*M_PI;
        sc.phi   = (output_width-x)/output_width*(2.0f*M_PI);
        
        s_coor3 coor3 = sphericalCoor2coor3(sc);
        
        world2cam(coor2, coor3);
        
        float yy = coor2.x;
        float xx = coor2.y;

        ocam_model.mappingConversionY(y,x) = yy;
        ocam_model.mappingConversionX(y,x) = xx;
        
      }
    }
  }
  
}
void Shc::init_lib_noise() {
  
  uniform_real_distribution<float> distribution(-1,1);
  bool bool_unsufficient = false;  
  int n_noise = n_bands;

  switch (noise_type) {
    case ZERO:
      noise_amplitudes.resize(n_bands);
      noise_amplitudes.setZero();
      break;
    case CONSTANT:
      noise_amplitudes.resize(n_bands);
      noise_amplitudes.setOnes();
      break;
    case NATURAL:
      noise_amplitudes.resize(n_bands);
      for (int l=0; l<n_bands; l++) {
        // formula obtained by fitting a function to a total of 75 amplitude spectra collected in indoor and outdoor environments
        noise_amplitudes(l) = 1.0f / (0.9983 + 0.9673 * l + 0.0013 * l*l); 
      }
      break;
    case CUSTOM:
      n_noise = noise_amplitudes.size();
      if (n_bands > n_noise) {
        bool_unsufficient = true;
      }
      break;    
  }
  
  noise_surf.resize(noise_samples);
  print_progress_start();
  for (int i=0; i<noise_samples; i++) {
    
    Coef coef(n_bands*n_bands);
    coef.setZero();
    
    if (noise_amplitudes.sum() > 0) {
      for (int l=0; l<n_bands; l++) {
        
        if (l >= n_noise) {
          break;
        }
        
        VectorReal vm(2*l+1);
        
        vm.setZero();
        for (int m=0; m<2*l+1; m++) {
          vm(m) = distribution(rand_generator);
        }
        
        float norm = vm.norm();
        if (norm == 0) {
          norm = 1;
        }
        
        float f = noise_amplitudes(l) / norm;
        for (int m=0; m<2*l+1; m++) {
          vm(m) *= f;
        }
        coef.segment(l*l, 2*l+1) = vm;    
      }
    }
    
    Shpm shpm = coef2shpm(coef);
    noise_surf[i] = intern_shpm2surf(shpm);
    
    float min = noise_surf[i].minCoeff();
    noise_surf[i] = noise_surf[i].array() - min;
    float mean = noise_surf[i].mean();
    if (mean == 0) {
      mean = 1;
    }
    noise_surf[i] = noise_surf[i].array() / mean;
    
    print_progress_update((float)(i) / (float)(noise_samples));
    
  }
  print_progress_stop();
  
  if (bool_unsufficient == true) {
    print_warning("init_lib_noise","the amplitude spectrum used to calculate the noise has less bands than the number of bands used.");
  }
  
  set_noise_enable(true);
  set_noise_mask(M_PI/2.0f);

}
bool Shc::check_init(bool is_init) {
  
  if (is_init == true) {
    
    if (bool_init == true) {
      print_warning("some_init_function", "Initialization already performed. No further calls to init functions possible!");
    }
    
    return !bool_init;
    
  } else {
  
    if (bool_init == false) {
      print_warning("some_none_init_function", "Shc object not initialized yet. Call init() first!");
    }
    
    return bool_init;
  
  }
    
}
bool Shc::check_init(Shpm& shpm) {
  
  if (shpm.initialized() == false) {
    print_warning("some_public_function", "A not initialized or empty Shpm instance has been passed to a function.");
    return false;
  }
  
  return true;
  
}
bool Shc::check_init(VecShpm& v_shpm) {
  
  for (uint i=0; i<(uint)v_shpm.size(); i++) {
    if (v_shpm[i].initialized() == false) {
      print_warning("some_public_function", "A not initialized or empty Shpm instance has been passed to a function.");
      return false;
    }
  }
  return true;
  
}
void Shc::init_rotations_none() {
  
  if (check_init(true) == false) {return;}
  
  s_rotation result;
  
  result.n_total_rotations = 0;
  result.l_max = 0;
  result.l_max_compass = 0;
  
  this->v_rotation.push_back(result);
  bool_init_rotations = true;
  
}
