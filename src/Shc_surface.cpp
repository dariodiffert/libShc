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

double Shc::calc_sh_P(int l, int m, double x) {   
  // evaluate an Associated Legendre Polynomial P(l,m,x) at x
  double pmm = 1;
  
  if(m>0) {
    double somx2 = sqrt((1-x)*(1+x));
    double fact = 1;
    for(int i=1; i<=m; i++) {
      pmm *= (-fact) * somx2;
      fact += 2;
    }
  }
  
  if(l==m) return pmm;
    
  double pmmp1 = x * (2*m+1) * pmm;
  
  if(l==m+1) return pmmp1;
  
  double pll = 0;
  
  for(int ll=m+2; ll<=l; ++ll) {
    pll = ( (2*ll-1)*x*pmmp1-(ll+m-1.0)*pmm ) / (ll-m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  
  return pll;
}
double Shc::calc_sh_K(int l, int m) { 
  
  double val; int am; double fac = 1;
  
  am = fabs(m);
  
  for (int f = l-am+1; f <= l+am; f++) {
    fac *= f;
  }
  
  val = sqrt((2*l+1)/(4*M_PI) / fac);
  
  return val;
} 
double Shc::calc_sh_Y(int l, int m, double theta, double phi) {
  
  double res;
  
  if (m > 0) {
    res = sqrt(2) * calc_sh_K(l,m) * cos( m*phi) * calc_sh_P(l, m,cos(theta));
  } else if (m < 0) {
    res = sqrt(2) * calc_sh_K(l,m) * sin(-m*phi) * calc_sh_P(l,-m,cos(theta));
  } else {
    res = calc_sh_K(l,0) * calc_sh_P(l,0,cos(theta));
  }
  
  return res;
  
}

vector<int> Shc::get_nearest_indices(s_sphericalCoor& sphericalCoor) {
  
  float theta = angular_normalization(sphericalCoor.theta);
  float phi   = angular_normalization(sphericalCoor.phi);
  
  if (theta > M_PI) {
    theta = 2*M_PI - theta;
    phi = angular_normalization(phi+M_PI);
  }
  
  int i1 = floor((theta + (lut_surf_quickRef.theta_stepSize/2.0)) / lut_surf_quickRef.theta_stepSize) - 1;
  int i2 = i1+1;
  
  if (i1 == -1) {
    i1 = 0;
    i2 = 0;
  }
  if (i2 >= lut_surf_quickRef.n_theta) {
    i1 = lut_surf_quickRef.n_theta-1;
    i2 = lut_surf_quickRef.n_theta-1;
  }
  
  int j1 = lut_surf_quickRef.theta_index[i1];
  int j2 = lut_surf_quickRef.theta_index[i2];
  
  int k11; int k12; int k21; int k22;
  k11 = floor((phi / (2.0*M_PI)) * lut_surf_quickRef.psi_count[i1]);
  if (k11 >= lut_surf_quickRef.psi_count[i1]-1) {
    k12 = 0;
  } else {
    k12 = k11+1;
  }
  k21 = floor((phi / (2.0*M_PI)) * lut_surf_quickRef.psi_count[i2]);
  if (k21 >= lut_surf_quickRef.psi_count[i2]-1) {
    k22 = 0;
  } else {
    k22 = k21+1;
  }
  
  vector<int> indices;

  indices.push_back(j1+k11);
  indices.push_back(j1+k12);
  indices.push_back(j2+k21);
  indices.push_back(j2+k22);
  
  return indices;
  
}
vector<float> Shc::get_nearest_weights(s_sphericalCoor& sphericalCoor, vector<int>& indices) {
  
  vector<float> result(4);
  
  float p, q;
  int i0, i1, i2, i3;
  float phi_min, phi_max;
  float theta_min, theta_max;
  
  i0 = indices[0]; i1 = indices[1]; i2 = indices[2]; i3 = indices[3];

  phi_min = min(lut_surf.unit[i0].sphericalCoor.phi, lut_surf.unit[i2].sphericalCoor.phi);
  phi_max = max(lut_surf.unit[i1].sphericalCoor.phi, lut_surf.unit[i3].sphericalCoor.phi);

  theta_min = lut_surf.unit[i0].sphericalCoor.theta;
  theta_max = lut_surf.unit[i2].sphericalCoor.theta;

  if (sphericalCoor.phi > phi_max) {sphericalCoor.phi = phi_max;}
  if (sphericalCoor.phi < phi_min) {sphericalCoor.phi = phi_min;}
  p = (sphericalCoor.phi - phi_min) / (phi_max - phi_min);

  if ((theta_max - theta_min) > 0) {
    q = (sphericalCoor.theta - theta_min) / (theta_max - theta_min);
  } else {
    q = 0;
  }

  result[0] = (1-p) * (1-q);
  result[1] = p * (1-q);
  result[2] = (1-p) * q;
  result[3] = p * q;
  
  return result;
      
}
float Shc::get_nearest_surf(Surf& surf, s_sphericalCoor& sphericalCoor) {
    
  vector<int> indices = get_nearest_indices(sphericalCoor);
  vector<float> weight = get_nearest_weights(sphericalCoor, indices);
    
  return Shc::get_nearest_surf(surf, indices, weight);
  
}
float Shc::get_nearest_surf(Surf& surf, vector<int>& indices, vector<float>& weight) {
  
  int n = indices.size();
  float weightSum = 0;
  for (int i=0; i<n; i++) {
    weightSum += weight[i];
  }
  
  float result = 0;
  for (int i=0; i<n; i++) {
    result += weight[i] * surf(indices[i]);
  }
  result /= weightSum;
  
  return result;
  
}
Shc::s_surf_quickRef Shc::create_surface_quickRef(s_surf& surface) {

  int n = surface.unit.size();
  
  s_surf_quickRef quickRef;
  quickRef.theta_stepSize = surface.unit[0].sphericalCoor.theta * 2;
  
  float last_theta = -1;
  int n_theta = 0;
  for (int i=0; i<n; i++) {
    if (last_theta != surface.unit[i].sphericalCoor.theta) {
      last_theta = surface.unit[i].sphericalCoor.theta;
      quickRef.theta.push_back(last_theta);
      quickRef.theta_index.push_back(i);
      n_theta++;
    }
  }
  
  quickRef.n_theta = n_theta;
  quickRef.psi_count.resize(n_theta);
  
  
  for (int i=0; i<quickRef.n_theta-1; i++) {
    quickRef.psi_count[i] = quickRef.theta_index[i+1]-quickRef.theta_index[i];
  }
  quickRef.psi_count[quickRef.n_theta-1] = n-quickRef.theta_index[quickRef.n_theta-1];
  
  return quickRef;
  
}
Shc::s_surf Shc::create_sphere(int n_points) {

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // n_points  = number of points equally distributed over the sphere (in the range theta_min to theta_max)
  //
  // adapted from: http://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int n_count;
  float theta; float phi;
  int M_theta; int M_phi; 
  float d_theta; float d_phi;
  float d;
  
  // area per point on sphere
  float areaPerAngle = 4 * M_PI / n_points;
  // assume that each area is a square, then d equals the altitude diameter of it
  d = sqrt(areaPerAngle);
  // total number of altitude steps for full sphere
  M_theta = round(M_PI/d);
  
  if (M_theta%2 == 0)
    M_theta ++;
  
  // altitude step size
  d_theta = M_PI / M_theta;
  // azimuth step size such that d_phi*d_theta=a
  d_phi   = areaPerAngle / d_theta;

  
  // go for it (calculate total number of entries)
  n_count = 0;
  for (int i=0; i<M_theta; i++) {
    
    theta = M_PI * (i+0.5)/M_theta;  
    
    M_phi = round(2*M_PI*sin(theta)/d_phi);
    if (M_phi%2 == 1)
      M_phi++;
    
    n_count += M_phi;
    
  }
  
  // allocate
  s_surf surface;
  surface.unit.resize(n_count);
  
  // go for it (calculate angles)
  n_count = 0;
  for (int i=0; i<M_theta; i++) {
    
    theta = M_PI * (i+0.5)/M_theta;
    
    M_phi = round(2*M_PI*sin(theta)/d_phi);
    if (M_phi%2 == 1)
      M_phi++;
    
    for (int j=0; j<M_phi; j++) {
      phi = 2*M_PI*j / M_phi;
      
      surface.unit[n_count].sphericalCoor.theta = theta;
      surface.unit[n_count].sphericalCoor.phi   = phi;
      surface.unit[n_count].weight              = 1;
      surface.unit[n_count].coor                = sphericalCoor2coor3(surface.unit[n_count].sphericalCoor);
      
      n_count++;
    }
  }
  
  float sum = 0;
  for (int i=0; i<n_count; i++) {
    sum += surface.unit[i].weight;
  }
  for (int i=0; i<n_count; i++) {
    surface.unit[i].weight = surface.unit[i].weight / sum * rescale_factor;
  }
  
  surface.n_angles = surface.unit.size();
  surface.n_sha    = l2sh(n_bands, FULL)*surface.n_angles;
  surface.l_max    = n_bands;
  surface.type     = SURF_UNIFORM_SPHERE;
  
  return surface;
  
}
Shc::s_surf Shc::create_sphere(int width, int height) {
  
  int n_count = width * height;
  
  // allocate
  s_surf surface;
  surface.unit.resize(n_count);
  
  // go for it (calculate angles)
  n_count = 0;
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
    
      float theta = y        * M_PI/ height;
      float phi   = x * 2.0f * M_PI/ width;
      
      surface.unit[n_count].sphericalCoor.theta = theta;
      surface.unit[n_count].sphericalCoor.phi   = phi;
      surface.unit[n_count].weight              = sin(theta);
      surface.unit[n_count].coor                = sphericalCoor2coor3(surface.unit[n_count].sphericalCoor);
      n_count++;
      
    }
  }
  
  float sum = 0;
  for (int i=0; i<n_count; i++) {
    sum += surface.unit[i].weight;
  }
  for (int i=0; i<n_count; i++) {
    surface.unit[i].weight = surface.unit[i].weight / sum * rescale_factor;
  }
  
  surface.n_angles = surface.unit.size();
  surface.n_sha    = l2sh(n_bands, FULL)*surface.n_angles;
  surface.l_max    = n_bands;
  surface.type     = SURF_UNIFORM_PANORAMA;
  
  return surface;
  
}
