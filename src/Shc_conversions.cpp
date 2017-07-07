//
// Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
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


Shc::s_sphericalCoor Shc::coor22sphericalCoor(s_coor2 coor2) {
  
  s_sphericalCoor sc;
  
  sc.phi   = angular_normalization(coor2.x * 2 * M_PI);
  sc.theta = angular_normalization(coor2.y * M_PI);
    
  return sc;
  
}
Shc::s_coor2 Shc::sphericalCoor2coor2(s_sphericalCoor sphericalCoor) {
  
  float t = angular_normalization(sphericalCoor.theta);
  float p = angular_normalization(sphericalCoor.phi);
  
  s_coor2 coor2;
  
  if (t <= M_PI) {
    coor2.y = t / M_PI;
    coor2.x = p / (2*M_PI);
  } else {
    coor2.y = (2.0*M_PI - t) / M_PI;
    coor2.x = (p+M_PI) / (2*M_PI);
  }
  
  return coor2;
  
}
Shc::s_coor3 Shc::sphericalCoor2coor3(s_sphericalCoor sphericalCoor) {
  
  s_coor3 coor3;
  
  coor3.x = sin(sphericalCoor.theta) * cos(sphericalCoor.phi);
  coor3.y = sin(sphericalCoor.theta) * sin(sphericalCoor.phi);
  coor3.z = cos(sphericalCoor.theta);
  
  return coor3;
  
}
Shc::s_sphericalCoor Shc::coor32sphericalCoor(s_coor3 coor3) {
  
  s_sphericalCoor sphericalCoor;
  
  float r = sqrt(coor3.x*coor3.x+coor3.y*coor3.y+coor3.z*coor3.z);
  
  if (r > 0) {
    sphericalCoor.theta = acos(coor3.z/r);
    sphericalCoor.phi   = atan2(coor3.y,coor3.x);
  } else {
    sphericalCoor.theta = M_PI/2.0f;
    sphericalCoor.phi   = 0;
  }

  return sphericalCoor;
  
}

void Shc::world2cam(s_coor2& coor2, s_coor3& coor3) {
  // coor3->coor2
  // adapted from Scaramuzza (https://sites.google.com/site/scarabotix/ocamcalib-toolbox)
  
  float inv_x=1.0f, inv_y=1.0f, inv_z=1.0f;
  if (ocam_model.invert_x_axis) {inv_x = -1.0f;}
  if (ocam_model.invert_y_axis) {inv_y = -1.0f;}
  if (ocam_model.invert_z_axis) {inv_z = -1.0f;}
  
  float norm  = sqrt(coor3.x*coor3.x + coor3.y*coor3.y + coor3.z*coor3.z);
  float theta = atan(inv_z * coor3.z/norm);
  float t, t_i, rho, x, y, invnorm;
  
  if (norm != 0) {
    invnorm = 1.0f/norm;
    t   = theta;
    rho = ocam_model.invpol[0];
    t_i = 1;

    for (int i=1; i<ocam_model.length_invpol; i++) {
      t_i *= t;
      rho += t_i*ocam_model.invpol[i];
    }

    x = inv_x * coor3.x*invnorm*rho;
    y = inv_y * coor3.y*invnorm*rho;
  
    coor2.x = x*ocam_model.c + y*ocam_model.d + ocam_model.xc;
    coor2.y = x*ocam_model.e + y              + ocam_model.yc;
  } else {
    coor2.x = ocam_model.xc;
    coor2.y = ocam_model.yc;
  }
  
}
void Shc::cam2world(s_coor2& coor2, s_coor3& coor3) {
  // coor2->coor3
  // adapted from Scaramuzza (https://sites.google.com/site/scarabotix/ocamcalib-toolbox)
  
  float inv_x=1.0f, inv_y=1.0f, inv_z=1.0f;
  if (ocam_model.invert_x_axis) {inv_x = -1.0f;}
  if (ocam_model.invert_y_axis) {inv_y = -1.0f;}
  if (ocam_model.invert_z_axis) {inv_z = -1.0f;}
  
  float invdet  = 1.0f/(ocam_model.c-ocam_model.d*ocam_model.e); // 1/det(A), where A = [c,d;e,1] as in the Matlab file

  float xp = invdet*(               (coor2.x - ocam_model.xc) - ocam_model.d*(coor2.y - ocam_model.yc) );
  float yp = invdet*( -ocam_model.e*(coor2.x - ocam_model.xc) + ocam_model.c*(coor2.y - ocam_model.yc) );

  float r   = sqrt(  xp*xp + yp*yp ); //distance [pixels] of  the point from the image center
  float zp  = ocam_model.pol[0];
  float r_i = 1;

  for (int i=1; i<ocam_model.length_pol; i++) {
    r_i *= r;
    zp  += r_i*ocam_model.pol[i];
  }

  //normalize to unit norm
  float invnorm = 1.0f/sqrt( xp*xp + yp*yp + zp*zp );

  coor3.x = inv_x * invnorm*xp;
  coor3.y = inv_y * invnorm*yp; 
  coor3.z = inv_z * invnorm*zp;

}

int Shc::getClosestIndex(float angle, int index_max, VectorReal& angles, vector<int>& fastIndices) {
  
  float d_cur;
  float d_min = 1e12;
  int minIndex = 0;
  
  for (int i=0; i<index_max; i++) {
    d_cur = angular_difference(angle, angles[fastIndices[i]]);
    if (d_cur < d_min) {
      d_min = d_cur;
      minIndex = i;
    }
  }
  
  return minIndex;
  
}
void Shc::xyz2rotationIndices(Xyz xyz, int index_rot_par, s_rotation_indices& indices) {
  
  indices.i = getClosestIndex(xyz.z, v_rotation[index_rot_par].n_fast_indices_z, v_rotation[index_rot_par].psi, v_rotation[index_rot_par].fast_indices_z);
  indices.j = getClosestIndex(xyz.y, v_rotation[index_rot_par].n_fast_indices_y, v_rotation[index_rot_par].psi, v_rotation[index_rot_par].fast_indices_y);
  indices.k = getClosestIndex(xyz.x, v_rotation[index_rot_par].n_fast_indices_x, v_rotation[index_rot_par].psi, v_rotation[index_rot_par].fast_indices_x);
  
  return;
  
}
Xyz Shc::rotationIndices2xyz(int index_rot_par, s_rotation_indices& indices) {
  
  Xyz xyz;
  
  xyz.z = v_rotation[index_rot_par].psi[v_rotation[index_rot_par].fast_indices_z[indices.i]];
  xyz.y = v_rotation[index_rot_par].psi[v_rotation[index_rot_par].fast_indices_y[indices.j]];
  xyz.x = v_rotation[index_rot_par].psi[v_rotation[index_rot_par].fast_indices_x[indices.k]];
  
  return xyz;
    
}

Shpm Shc::intern_surf2shpm_hemi_FFT(Surf& surf, e_contin contin, int l_max) {
  
  Shpm result;
  result.coef.resize(l2sh(l_max, FULL));
  result.coef.setZero();
  result.contin = FULL;
  result.l_max = l_max;
  
  int l1=0, l2=0;
 
  switch (contin) {
    case HEMI_M:
      l1 = 0; l2 = 1;
      break;
    case HEMI_MN:
      l1 = 1; l2 = 0;
      break;
    case HEMI_RM:
      l1 = 0; l2 = 0;
      break;
    case HEMI_RMN:
      l1 = 1; l2 = 1;  
      break;
    default:
      print_warning("intern_surf2shpm_hemi", "this should not be reachable. BUG!?");
  }
  
  int n_mid = (lut_surf_quickRef.n_theta-1)/2;
  for (int t=0; t<=n_mid; t++) {
    int x = lut_FFT.theta_index[t];
    int s = lut_surf_quickRef.theta_index[t];
    int n = lut_FFT.store[t].size();
    for (int i=0; i<n; i++) { // since the size of surf is knwon, this could be shortened !!!
      lut_FFT.store[t][i].r = 0;
      lut_FFT.store[t][i].i = 0;
    }
    lut_FFT.cfg_FFT[x].compute(surf, s, lut_FFT.store[t]);
  }
  
  int ll;
  for (int l=0; l<l_max; l++) {
    ll = l%2;
    if ((ll == l1) || (ll == l2) || l == 0) {
    
      int i = sh2index(l,0);
      for (int t=0; t<n_mid; t++) {
        result.coef(i) += 2.0 * lut_FFT.val_FFT[i*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][0].r;
      }
      
      for (int m=1; m<=l; m++) {
        int i1 = sh2index(l, m);
        int i2 = sh2index(l,-m);
        for (int t=0; t<n_mid; t++) {
          result.coef(i1) += 2.0 * lut_FFT.val_FFT[i1*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][m].r;
          result.coef(i2) += 2.0 * lut_FFT.val_FFT[i2*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][m].i;
        }      
      }
      
    }
  }
  
  return result;
  
}
Shpm Shc::intern_surf2shpm_hemi(Surf& surf, e_contin contin, int l_max) {
  
  Shpm result;
  result.coef.resize(l2sh(l_max, FULL));
  result.coef.setZero();
  result.contin = FULL;
  result.l_max = l_max;
  
  int m1=0; int l1=0;
  int m2=0; int l2=0;
 
  switch (contin) {
    case HEMI_M:
      m1 = 0; l1 = 0;
      m2 = 1; l2 = 1;
      break;
    case HEMI_MN:
      m1 = 0; l1 = 1;
      m2 = 1; l2 = 0;
      break;
    case HEMI_RM:
      m1 = 0; l1 = 0;
      m2 = 1; l2 = 0;
      break;
    case HEMI_RMN:
      m1 = 0; l1 = 1;
      m2 = 1; l2 = 1;  
      break;
    default:
      print_warning("intern_surf2shpm_hemi", "this should not be reachable. BUG!?");
  }
  
  int ll, mm;
  int count = 0;
  for (int l=0; l<l_max; l++) { 
    ll = l%2;
    for (int m=-l; m<=l; m++) {
      mm = (m%2+2)%2;
      if ((ll == l1 && mm == m1) || (ll == l2 && mm == m2) || l == 0) {
        int i = sha2index(lut_surf,l,m,0);
        result.coef[count] = surf.dot(lut_surf.sh_weighted_hemi.segment(i, lut_surf.n_angles));
      }
      count++;
    }
  }  
  
  return result;
  
}
Shpm Shc::intern_surf2shpm_full_FFT(Surf& surf, int l_max) {
    
  Shpm result;
  result.coef.resize(l2sh(l_max, FULL));
  result.coef.setZero();
  result.contin = FULL;
  result.l_max = l_max; 
  
  for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
    int x = lut_FFT.theta_index[t];
    int s = lut_surf_quickRef.theta_index[t];
    int n = lut_FFT.store[t].size();
    for (int i=0; i<n; i++) { // since the size of surf is knwon, this could be shortened !!!
      lut_FFT.store[t][i].r = 0;
      lut_FFT.store[t][i].i = 0;
    }
    lut_FFT.cfg_FFT[x].compute(surf, s, lut_FFT.store[t]);
  }
  
  for (int l=0; l<l_max; l++) {
    
    int i = sh2index(l,0);
    for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
      result.coef(i) += lut_FFT.val_FFT[i*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][0].r;
    }
    
    for (int m=1; m<=l; m++) {
      int i1 = sh2index(l, m);
      int i2 = sh2index(l,-m);
      for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
        result.coef(i1) += lut_FFT.val_FFT[i1*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][m].r;
        result.coef(i2) += lut_FFT.val_FFT[i2*lut_surf_quickRef.n_theta+t] * lut_FFT.store[t][m].i;
      }      
    }
  }
  
  return result;
  
} 
Shpm Shc::intern_surf2shpm_full(Surf& surf, int l_max) {
  
  Shpm result;
  result.coef.resize(l2sh(l_max, FULL));
  result.coef.setZero();
  result.contin = FULL;
  result.l_max = l_max;
  
  int count = 0;
  for (int l=0; l<l_max; l++) {    
    for (int m=-l; m<=l; m++) {
      int i = sha2index(lut_surf,l,m,0);
      result.coef[count] = surf.dot(lut_surf.sh_weighted_full.segment(i, lut_surf.n_angles));
      count++;
    }
  }

  return result;
  
}
Shpm Shc::intern_surf2shpm(Surf& surf, e_contin contin, int l_max, bool allow_noise, e_contin contin_result) {
   
  Shpm shpm;
  Surf surf_mod(lut_surf.n_angles);
  
  if (allow_noise == true && noise_enable == true && noise_samples > 0 && contin == FULL) {
  
    uniform_int_distribution<int> distribution(0,noise_samples-1);
    int noise_i = distribution(rand_generator);
    
    Surf MI = (1.0f-noise_mask.array()) * surf.array();
    Surf MN =       noise_mask.array()  * noise_surf[noise_i].array();
    
    float mi = MI.sum();
    float mn = MN.sum();
    if (lut_surf.n_angles-noise_mask_sum != 0) { 
      mi /= (lut_surf.n_angles-noise_mask_sum);
    }
    if (noise_mask_sum != 0) {
      mn /= noise_mask_sum;
    }
    
    float si=0;
    for (int i=0; i<lut_surf.n_angles; i++) {
      si += ((MI[i]-mi)*(MI[i]-mi)) * (1.0f-noise_mask[i]);
    }
    si /= (lut_surf.n_angles-noise_mask_sum);

    float sn=0;
    for (int i=0; i<lut_surf.n_angles; i++) {
      sn += ((MN[i]-mn)*(MN[i]-mn)) * noise_mask[i];
    }
    sn /= noise_mask_sum;

    if (mn == 0 || sn == 0) {
      for (int i=0; i<lut_surf.n_angles; i++) {
        surf_mod[i] = (noise_mask[i] * mn * noise_amplifier + (1.0f-noise_mask[i]) * surf[i]);
      }
    } else {
      float a = sqrt(si/sn)*mi/mn;
      float b = mi*(1.0f-sqrt(si/sn));
      for (int i=0; i<lut_surf.n_angles; i++) {
        // float val2 = ((noise_surf[noise_i][i] - mn) * sqrt(si/sn) + mn) * mi/mn; //substituting a and b simplifies calculation
        float val = a * noise_surf[noise_i][i] + b;
        surf_mod[i] = (noise_mask[i] * val * noise_amplifier + (1.0f-noise_mask[i]) * surf[i]);
      }
    }
     
  } else {
    surf_mod = surf;
  }
  
  switch (contin) {
    case FULL:
      if (bool_FFT) {
        shpm = intern_surf2shpm_full_FFT(surf_mod, l_max); 
      } else {
        shpm = intern_surf2shpm_full(surf_mod, l_max);
      }
      break;
    default:
      if (bool_FFT) {
        shpm = intern_surf2shpm_hemi_FFT(surf_mod, contin, l_max); 
      } else {
        shpm = intern_surf2shpm_hemi(surf_mod, contin, l_max);
      }
      break;    
  }
  
  if (shpm.contin != contin_result) {
    shpm = contin_convert(shpm, contin_result);
  }
  
  return shpm;
  
}

Surf Shc::intern_shpm2output(Shpm& shpm) {
  
  Coef coef = shpm2coef(shpm);
  int l_max = shpm.l_max;
  
  VectorReal ifft(output_width);
  Surf result(output_width*output_height);
  result.setZero();
  
  for (int h=0; h<output_height; h++) {
    for (int l=0; l<l_max; l++) {
      
      int w=0;
      int i = sh2index(l, 0);
      lut_IFFT_output.in[w].r = coef(i) * lut_IFFT_output.val[i * output_height + h];
      lut_IFFT_output.in[w].i = 0;
      w++;
      
      for (int m=1; m<=l; m++) {
        
        if (w >= output_width) {break;}
        
        int i1 = sh2index(l, m);
        int i2 = sh2index(l,-m);
        lut_IFFT_output.in[w].r = coef(i1) * lut_IFFT_output.val[i1 * output_height + h];
        lut_IFFT_output.in[w].i = coef(i2) * lut_IFFT_output.val[i2 * output_height + h];
        w++;

      }
      
      for (int wi=w; wi<output_width; wi++) {
        lut_IFFT_output.in[wi].r = 0;
        lut_IFFT_output.in[wi].i = 0;
      }
      
      lut_IFFT_output.cfg.compute(lut_IFFT_output.in, ifft, 0);
      result.segment(h*output_width, output_width) += ifft.reverse();
      
    }
  }
  
  return result;

}
Surf Shc::intern_shpm2surf(Shpm& shpm) {
  
  Coef coef = shpm2coef(shpm);
  int l_max = shpm.l_max;
  
  Surf result(lut_surf.n_angles);
  result.setZero();
  
  if (bool_sphere_mode == true && bool_FFT == true) {
    int p_count = 0;
    for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
      
      int p_max = lut_FFT.store[t].size();
      int p_cur = lut_surf_quickRef.psi_count[t];
      VectorReal ifft(p_max);

      for (int l=0; l<l_max; l++) {
        
        if (contin_skip(shpm.contin, l)) {continue;}
        
        int i = sh2index(l, 0);
        lut_FFT.store[t][0].r = coef(i) * lut_FFT.val_IFFT[i * lut_surf_quickRef.n_theta + t];
        lut_FFT.store[t][0].i = 0;
        
        for (int p=1; p<p_max; p++) {
          if (p <= l) {
            int i1 = sh2index(l, p);
            int i2 = sh2index(l,-p);
            lut_FFT.store[t][p].r = coef(i1) * lut_FFT.val_IFFT[i1 * lut_surf_quickRef.n_theta + t];
            lut_FFT.store[t][p].i = coef(i2) * lut_FFT.val_IFFT[i2 * lut_surf_quickRef.n_theta + t];
          } else {
            lut_FFT.store[t][p].r = 0;
            lut_FFT.store[t][p].i = 0;
          }
        }
        
        int x = lut_FFT.theta_index[t];
        lut_FFT.cfg_IFFT[x].compute(lut_FFT.store[t], ifft, 0);
        result.segment(p_count, p_cur) += ifft.segment(0, p_cur);
        
      }
      p_count += p_cur;
    }
  } else {
    MatrixReal matrix = shpm2matrix(shpm);
    result = matrix2surf(matrix);
    print_warning("intern_shpm2surf", "Inverse Fourier transform on non-spherical (custom) surface; output is used as surface instead. This can happen if noise is used together with non-spherical surface, in this case use sufficiently large output dimensions!");
  }
  
  return result;
  
}
Surf Shc::coef2surf(int index) {
  
  Surf result(lut_surf.n_angles);
  result.setZero();
  
  int l, m;
  index2sh(index, l, m);
  
  int p_count = 0;
  for (int t=0; t<lut_surf_quickRef.n_theta; t++) {
    
    int p_max = lut_FFT.store[t].size();
    int p_cur = lut_surf_quickRef.psi_count[t];
    VectorReal ifft(p_max);

    int i = sh2index(l, m);
    
    for (int p=0; p<p_max; p++) {
      lut_FFT.store[t][p].r = 0;
      lut_FFT.store[t][p].i = 0;
    }
    
    int m_abs = abs(m);
    if (m == 0) {
      lut_FFT.store[t][m_abs].r = lut_FFT.val_IFFT[i * lut_surf_quickRef.n_theta + t];
      lut_FFT.store[t][m_abs].i = 0;
    } else if (m > 0) {
      lut_FFT.store[t][m_abs].r = lut_FFT.val_IFFT[i * lut_surf_quickRef.n_theta + t];
      lut_FFT.store[t][m_abs].i = 0;
    } else if (m < 0) {
      lut_FFT.store[t][m_abs].r = 0;
      lut_FFT.store[t][m_abs].i = lut_FFT.val_IFFT[i * lut_surf_quickRef.n_theta + t];
    }
    
    int x = lut_FFT.theta_index[t];
    lut_FFT.cfg_IFFT[x].compute(lut_FFT.store[t], ifft, 0);
    result.segment(p_count, p_cur) += ifft.segment(0, p_cur);

    p_count += p_cur;
  }
  
  return result;
  
}

CoefBandwise Shc::shpm2bandwise(Shpm& shpm, int l_max) {
  
  vector<VectorReal> result(n_bands);
  int count = 0;
  
  if (l_max < 0)
    l_max = shpm.l_max;
  
  if (trifold(shpm.contin) == FULL) {
    
    // spherical case
    for (int l=0; l<n_bands; l++) {
      result[l].resize(2*l+1);
      if (l < l_max) {
        for (int m=-l; m<=l; m++) {
          result[l](m+l) = shpm.coef[count];       
          count++;
        }
      } else {
        result[l].setZero();
      }
    }
    
  } else {   
    
    // hemispherical case
    for (int l=0; l<n_bands; l++) {
      result[l].resize(2*l+1);
      if ((l >= l_max) || contin_skip(shpm.contin, l)) {
        result[l].setZero();
      } else {
        for (int m=-l; m<=l; m++) {
          result[l](m+l) = shpm.coef[count];
          count++;
        }
      }
      
    }
    
  }
  
  return result;
  
}

Shpm Shc::full2hemi(Shpm& shpm, e_contin contin) {
    
  Coef hemi(l2sh(shpm.l_max, contin));
  hemi.setZero();
  hemi(0) = shpm.coef(0);
  
  int count_full = 1;
  int count_hemi = 1;
  for (int l=1; l<shpm.l_max; l++) {
    
    int n = 2*l+1;
    if ((contin == HEMI_RM && l%2 == 0) || (contin == HEMI_RMN && l%2 == 1)) {
      hemi.segment(count_hemi, n) = shpm.coef.segment(count_full, n);
      count_hemi += n;
    }
    count_full += n;
    
  }
  
  Shpm result;
  result.coef   = hemi;
  result.contin = contin;
  result.l_max  = shpm.l_max;
  
  return result;
  
}
Shpm Shc::hemi2full(Shpm& shpm) {
  
  Coef full(l2sh(shpm.l_max, FULL));
  full.setZero();
  full(0) = shpm.coef(0);
  
  int count_full = 1;
  int count_hemi = 1;
  for (int l=1; l<shpm.l_max; l++) {
    
    int n = 2*l+1;
    if ((shpm.contin == HEMI_RM && l%2 == 0) || (shpm.contin == HEMI_RMN && l%2 == 1)) {
      full.segment(count_full, n) = shpm.coef.segment(count_hemi, n);
      count_hemi += n;
    }
    count_full += n;
    
  }
  
  Shpm result;
  result.coef   = full;
  result.contin = FULL;
  result.l_max  = shpm.l_max;
  
  return result;
  
}

MatrixReal Shc::intern_output2matrix(Surf& surf) {
  
  if (!check_init()) {return MatrixReal();}
  
  MatrixReal out;
  out.resize(output_height, output_width);
  
  int count = 0;
  for (int y=0; y<output_height; y++) {
    for (int x=0; x<output_width; x++) {
      
      out(y,x) = surf[count];
      count++;
      
    }
  }
  
  return out;
    
}

Shpm Shc::contin_convert(Shpm& shpm, e_contin contin, int l_max) {
  
  if (shpm.coef.size() != l2sh(shpm.l_max, shpm.contin)) {
    print_warning("contin_convert", "Input vector has invalid size! BUG");
  }
  
  Shpm result;
  
  int l_max_old = shpm.l_max;
  if (l_max == -1 || l_max > shpm.l_max) {
    l_max = shpm.l_max;
  }
  shpm.l_max = l_max;
    
  if (trifold(shpm.contin) == trifold(contin)) {
    result.contin = contin;
    result.l_max  = shpm.l_max;
    result.coef   = shpm.coef.segment(0, l2sh(shpm.l_max, contin));
  } else if (trifold(shpm.contin) == FULL) {
    result = full2hemi(shpm, contin);
  } else if (trifold(shpm.contin) != FULL && trifold(contin) != FULL) {
    result.contin = contin;
    result.l_max  = shpm.l_max;
    result.coef   = VectorReal(l2sh(result.l_max, result.contin)); result.coef.setZero();
  } else {
    result = hemi2full(shpm);
  }
  
  shpm.l_max = l_max_old;
  
  return result;

}
