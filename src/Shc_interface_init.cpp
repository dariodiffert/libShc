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


void Shc::init_tolerances(float translation, float cg) {
  
  if (!check_init(true)) {return;}
  
  if (translation < 0 || cg < 0) {
    print_warning("init_tolerances", "Negative tolerance values are not valid.");
    return;
  }
  
  tolerances.translation = translation;
  tolerances.cg          = cg;
  bool_init_tolerances = true;
  
}
void Shc::init_files(string file_surface, string file_translation, bool force_overwrite) {
  
  if (!check_init(true)) {return;}
  
  files.surface         = file_surface;
  files.translation     = file_translation;
  files.force_overwrite = force_overwrite;
  bool_init_files = true;
  
}
void Shc::init_rotations_sphere(float resolution, float range, int l_max_compass, int l_max) {
  
  if (!check_init(true)) {return;}

  float tilt_max = range / 2.0f;
  
  if (tilt_max > M_PI+1e-4 || tilt_max < 0) {
    print_warning("init_rotations_sphere", "range has to be in (0,2*pi].");
    return;
  }
  
  s_rotation result;
    
  result.psi = stepspace(resolution, 2*M_PI-resolution, false);
  result.theta.resize(2);
  result.theta[0] = -M_PI/2;
  result.theta[1] =  M_PI/2;
  result.phi.resize(2);
  result.phi[0] = -M_PI/2;
  result.phi[1] =  M_PI/2;
  VectorReal angles = stepspace(resolution, tilt_max, true); //linspace(-tilt_max, tilt_max, tilt_steps);
  create_fast_indices(angles, resolution, result);   
  
  result.n_total_rotations = result.fast_indices_x.size() * result.fast_indices_y.size() * result.fast_indices_z.size();
     
  result.l_max         = l_max;
  result.l_max_compass = l_max_compass;
  result.resolution    = resolution;
  
  this->v_rotation.push_back(result);
  bool_init_rotations = true;
  
}
void Shc::init_rotations_cone(float resolution, float range, int l_max_compass, int l_max) {
  
  if (!check_init(true)) {return;}
  
  float cone = range / 2.0f;
  
  if (cone > M_PI || cone <= 0) {
    print_warning("init_rotations_cone", "range has to be in (0,pi].");
    return;
  }
  
  s_rotation result; 
  
  result.psi = stepspace(resolution, cone, true);
  result.theta.resize(2);
  result.theta[0] = -M_PI/2;
  result.theta[1] =  M_PI/2;
  result.phi.resize(2);
  result.phi[0] = -M_PI/2;
  result.phi[1] =  M_PI/2;
  create_fast_indices(result.psi, resolution, result);

  result.n_total_rotations = result.fast_indices_x.size() * result.fast_indices_y.size() * result.fast_indices_z.size();
  
  result.l_max         = l_max;
  result.l_max_compass = l_max_compass;
  result.resolution    = resolution;
  
  this->v_rotation.push_back(result);
  bool_init_rotations = true;
  
}
void Shc::init_rotations_custom(VectorReal rotation_angles, VectorInt indices_x, VectorInt indices_y, VectorInt indices_z, int l_max_compass, int l_max) {
  
  if (!check_init(true)) {return;}
  
  int n = rotation_angles.size();
  if (n <= 0) {
    print_warning("init_rotations_custom", "No rotation angles passed.");
    return;
  }
  if (indices_x.size() == 0 || indices_y.size() == 0 || indices_z.size() == 0) {
    print_warning("init_rotations_custom", "At least one index has to be passed for each axis.");
    return;
  }
    
  s_rotation result; 
  
  result.psi = rotation_angles;
  result.theta.resize(2);
  result.theta[0] = -M_PI/2;
  result.theta[1] =  M_PI/2;
  result.phi.resize(2);
  result.phi[0] = -M_PI/2;
  result.phi[1] =  M_PI/2;
  
  result.n_fast_indices_z = indices_z.size();
  result.n_fast_indices_y = indices_y.size();
  result.n_fast_indices_x = indices_x.size();
  
  if (result.n_fast_indices_z == 0 || result.n_fast_indices_y == 0 || result.n_fast_indices_x == 0) {
    print_warning("init_rotations_custom", "All indices vectors have to contain at least one entry.");
    return;
  }
  
  for (int i=0; i<result.n_fast_indices_x; i++) {
    if (indices_x[i] < 0 || indices_x[i] >= n) {print_warning("init_rotations_custom", "indices_x: Invalid indices"); return;}
    result.fast_indices_x.push_back(indices_x[i]);
  }
  for (int i=0; i<result.n_fast_indices_y; i++) {
    if (indices_y[i] < 0 || indices_y[i] >= n) {print_warning("init_rotations_custom", "indices_y: Invalid indices"); return;}
    result.fast_indices_y.push_back(indices_y[i]);
  }
  for (int i=0; i<result.n_fast_indices_z; i++) {
    if (indices_z[i] < 0 || indices_z[i] >= n) {print_warning("init_rotations_custom", "indices_z: Invalid indices"); return;}
    result.fast_indices_z.push_back(indices_z[i]);
  }
    
  result.n_total_rotations = result.fast_indices_x.size() * result.fast_indices_y.size() * result.fast_indices_z.size();
  
  result.l_max         = l_max;
  result.l_max_compass = l_max_compass;
  result.resolution    = -1;
  
  this->v_rotation.push_back(result);
  bool_init_rotations = true;
  
}
void Shc::init_surface(int n_points, bool FFT) {
  
  if (!check_init(true)) {return;}
  
  if (bool_init_surface == true) {
    print_warning("init_surface", "Calling init_surface() multiple times is not allowed!");
    return;
  }
  
  if (n_points < 1) {
    print_warning("init_surface", "The number of sample points has to be greater zero.");
    return;
  }
    
  this->bool_sphere_mode = true;
  this->bool_FFT         = FFT;
  this->lut_surf         = create_sphere(n_points);
  
  bool_init_surface = true;

}
void Shc::init_surface(int width, int height, bool FFT) {
  
  if (!check_init(true)) {return;}
  
  if (bool_init_surface == true) {
    print_warning("init_surface", "Calling init_surface() multiple times is not allowed!");
    return;
  }
  
  if (width < 1 || height < 1) {
    print_warning("init_surface", "width and height have to be values greater zero.");
    return;
  }
  
  if (width % 2 == 1 && FFT == true) {
    print_warning("init_surface", "FFT can only be used with an even width.");
    return;
  }
    
  this->bool_sphere_mode = true;
  this->bool_FFT         = FFT;
  this->lut_surf         = create_sphere(width, height);
  
  bool_init_surface = true;

}
void Shc::init_surface(VectorReal theta, VectorReal phi, VectorReal weight) {
  
  if (!check_init(true)) {return;}
  
  if (bool_init_surface == true) {
    print_warning("init_surface", "Calling init_surface() multiple times is not allowed!");
    return;
  }
  
  s_surf surface;
  
  int n = phi.size();
  if (n != (int)theta.size() || n != (int)weight.size()) {
    print_warning("init_surface", "All passed vectors have to be of the same size.");
    return;
  }
  if (n == 0) {
    print_warning("init_surface", "The number of sample points has to be greater zero.");
    return;
  }
  
  surface.unit.resize(n);
  
  for (int i=0; i<n; i++) {
    surface.unit[i].sphericalCoor.theta = theta[i];
    surface.unit[i].sphericalCoor.phi   = phi[i];
    surface.unit[i].weight              = weight[i];
    surface.unit[i].coor                = sphericalCoor2coor3(surface.unit[i].sphericalCoor);
  }
  
  float sum = 0;
  for (int i=0; i<n; i++) {
    sum += surface.unit[i].weight;
  }
  for (int i=0; i<n; i++) {
    surface.unit[i].weight = surface.unit[i].weight / sum * rescale_factor;
  }
  
  surface.type = SURF_CUSTOM;
  
  this->bool_sphere_mode = false;
  this->bool_FFT         = false;
  this->lut_surf         = surface;
  bool_init_surface      = true;
  
}
void Shc::init_surface(string file_surface) {
  
  if (!check_init(true)) {return;}
  
  if (bool_init_surface == true) {
    print_warning("init_surface", "Calling init_surface() multiple times is not allowed!");
    return;
  }
  
  vector <vector <string> > data; 
  ifstream file(file_surface);

  if (file.is_open() == false) {
    print_warning("init_surface", "File not found!");
    return;
  }
  
  while (file) {
    string line;
    if (!getline(file, line))
      break;
    
    stringstream ss(line);
    vector <string> record;

    while (ss) {
      string line;
      if (!getline(ss, line, ','))
        break;
      record.push_back(line);
    }

    data.push_back(record);
  }
  
  int n = data.size();
  
  if (n <= 0) {
    print_warning("init_surface", "No data points found in file!");
    return;
  }
  
  VectorReal theta(n);
  VectorReal phi(n);
  VectorReal weight(n);
  
  for (int i=0; i<n; i++) {
    theta[i]  = atof(data[i][0].c_str());
    phi[i]    = atof(data[i][1].c_str());
    weight[i] = atof(data[i][2].c_str());
  }
  
  init_surface(theta, phi, weight);
  
}
void Shc::init_output(int width, int height) {
  
  if (!check_init(true)) {return;}
  
  if (width <= 0 || height <= 0) {
    print_warning("init_output", "Invalid dimensions, width and height need to be at least be set to 2.");
    return;
  }
  if (width%2==1) {
    print_warning("init_output", "Width has to be even.");
    width++;
  }
  if (height%2==1) {
    print_warning("init_output", "Height has to be even.");
    height++;
  }
  
  this->output_width  = width;
  this->output_height = height;
  
  bool_init_output = true;
  
}
void Shc::init_bands(int n_bands, int n_bands_CG) {
  
  if (!check_init(true)) {return;}
  
  if (n_bands <= 0) {
    print_warning("init_bands", "The number of used bands needs to be at least 1");
    return;
  }
  
  if (n_bands < n_bands_CG) {
    print_warning("init_bands", "The number of used bands for the Clebsch-Gordan coefficients can not exceed the number of maximal bands.");
    n_bands_CG = n_bands;
  }
  
  this->n_bands         = n_bands;
  this->n_bands_CG      = n_bands_CG;
  this->bool_init_bands = true;
  
}
void Shc::init_translations(int n_distances, e_trans_type trans_type, int l_max) {
  
  if (!check_init(true)) {return;}
  
  if (n_distances < 0) {
    print_warning("init_translations", "The number of translation steps (n_distances) has to be a positive number.");
    return;
  }
  
  VectorReal distances;
  if (n_distances > 1) {
    distances = linspace(0, 1, n_distances);
  }
  
  init_translations(distances, trans_type, l_max);
  
}
void Shc::init_translations(VectorReal distances, e_trans_type trans_type, int l_max) {
  
  if (!check_init(true)) {return;}
  
  int n = distances.size();
  
  for (int i=0; i<n; i++) {
    if (distances(i) < 0) {
      print_warning("init_translations", "Distances can only have positive values.");
      return;
    }
  }
  
  for (int i=0; i<n-1; i++) {
    if (distances(i) > distances(i+1)) {
      print_warning("init_translations", "Distances have to be in ascending order.");
      return;
    }
  }
  
  this->lut_translate.distances   = distances;
  this->lut_translate.n_distances = distances.size();
  this->lut_translate.l_max       = l_max;
  this->lut_translate.trans_type  = trans_type;
  
  sort(lut_translate.distances.data(), lut_translate.distances.data()+lut_translate.distances.size(), [](float lhs, float rhs){return rhs > lhs;});
    
  bool_init_translations = true;
  
}
void Shc::init_slices(int n_slices) {
  
  if (!check_init(true)) {return;}
  
  if (n_slices <= 0) {
    print_warning("init_slices", "Slices can only have positive values.");
    return;
  }

  VectorReal slices(n_slices);
  if (n_slices > 1) {
    VectorReal slices_temp = linspace(0,1,n_slices+1);
    for (int i=1; i<n_slices+1; i++) {
      slices(i-1) = slices_temp(i);
    }
  } else {
    slices(0) = 1;
  }
  
  init_slices(slices);
  
}
void Shc::init_slices(VectorReal slices) {
  
  if (!check_init(true)) {return;}
  
  int n_slices = slices.size();
  for (int i=0; i<n_slices; i++) {
    if (slices(i) <= 0) {
      print_warning("init_slices", "Slices can only have positive values.");
      return;
    }
  }
  
  for (int i=0; i<n_slices-1; i++) {
    if (slices(i) > slices(i+1)) {
      print_warning("init_slices", "Slices have to be in ascending order.");
      return;
    }
  }
  
  sort(slices.data(), slices.data()+slices.size(), [](float lhs, float rhs){return rhs > lhs;});
  
  this->lut_translate.slices   = slices;
  this->lut_translate.n_slices = n_slices;
  
  bool_init_slices = true;
  
}
void Shc::init_noise(int n_samples, VectorReal amplitude_spectrum) {
    
  if (!check_init(true)) {return;}
  
  if (n_samples < 0) {
    this->noise_samples = 0;
    print_warning("init_noise", "The number of samples cannot be negative");
    return;
  }
  
  for (uint i=0; i<(uint)amplitude_spectrum.size(); i++) {
    if (amplitude_spectrum[i] < 0) {
      print_warning("init_noise", "All entries of the amplitude spectrum have to be greater equal zero");
      return;
    }
  }
  
  this->noise_type       = CUSTOM;
  this->noise_amplitudes = amplitude_spectrum;
  this->noise_samples    = n_samples;
  this->bool_init_noise  = true;

}
void Shc::init_noise(int n_samples, e_noise noise) {
    
  if (!check_init(true)) {return;}
  
  if (n_samples < 0) {
    this->noise_samples = 0;
    print_warning("init_noise", "The number of samples cannot be negative");
    return;
  }
  
  if (noise == CUSTOM) {
    print_warning("init_noise", "To create noise based on a custom amplitude spectrum, use void Shc::init_noise(int n_samples, VectorReal amplitude_spectrum)");
    return;
  }
  
  this->noise_type       = noise;
  this->noise_amplitudes = VectorReal();
  this->noise_samples    = n_samples;
  this->bool_init_noise  = true;

}
void Shc::init_ocamcalib(string file, bool invert_x_axis, bool invert_y_axis, bool invert_z_axis) {
  // adapted from Scaramuzza (https://sites.google.com/site/scarabotix/ocamcalib-toolbox)
  
  if (!check_init(true)) {return;}
  
  if (file == "") {return;}

  const int CMV_MAX_BUF = 1024;
  FILE *f;
  char buf[CMV_MAX_BUF];
  int i;

  //Open file
  f=fopen(file.c_str(),"r");
  if (f == 0) {
    print_warning("init_ocamcalib", "File not found\n");                               
    return;
  }

  //Read polynomial coefficients
  (void)(fgets(buf,CMV_MAX_BUF,f)+1);
  (void)(fscanf(f,"\n")+1);
  (void)(fscanf(f,"%d", &ocam_model.length_pol)+1);
  for (i=0; i<ocam_model.length_pol; i++) {
    (void)(fscanf(f," %lf",&ocam_model.pol[i])+1);
  }

  //Read inverse polynomial coefficients
  (void)(fscanf(f,"\n")+1);
  (void)(fgets(buf,CMV_MAX_BUF,f)+1);
  (void)(fscanf(f,"\n")+1);
  (void)(fscanf(f,"%d", &ocam_model.length_invpol)+1);
  for (i=0; i<ocam_model.length_invpol; i++) {
    (void)(fscanf(f," %lf",&ocam_model.invpol[i])+1);
  }

  //Read center coordinates
  (void)(fscanf(f,"\n")+1);
  (void)(fgets(buf,CMV_MAX_BUF,f)+1);
  (void)(fscanf(f,"\n")+1);
  (void)(fscanf(f,"%lf %lf\n", &ocam_model.xc, &ocam_model.yc)+1);

  //Read affine coefficients
  (void)(fgets(buf,CMV_MAX_BUF,f)+1);
  (void)(fscanf(f,"\n")+1);
  (void)(fscanf(f,"%lf %lf %lf\n", &ocam_model.c, &ocam_model.d, &ocam_model.e)+1);

  //Read image size
  (void)(fgets(buf,CMV_MAX_BUF,f)+1);
  (void)(fscanf(f,"\n")+1);
  (void)(fscanf(f,"%d %d", &ocam_model.height, &ocam_model.width)+1);

  fclose(f);
  
  ocam_model.invert_x_axis = invert_x_axis;
  ocam_model.invert_y_axis = invert_y_axis;
  ocam_model.invert_z_axis = invert_z_axis;
  
  this->bool_init_ocamcalib  = true;
  
}

void Shc::init() { 
  
  #if !defined(EIGEN_VECTORIZE)
    print_warning("init","EIGEN_VECTORIZE is false. Make sure to activate Eigen vectorization instructions at compile time (see Eigen manual).");
  #endif
  
  // init default parameters if not specified otherwise
   
  print("-------- Initialization Started --------\n");
  
  if (bool_init_files == false) {
    print("init: Using standard values: init_files(\"\", \"\", \"\") (disabled).\n");
    init_files("", "", true);
  }
  
  if (bool_init_tolerances == false) {
    print("init: Using standard values: init_tolerances(0, 1e-12).\n");
    init_tolerances(0, 1e-12);
  }
  
  if (bool_init_bands == false) {
    print("init: Using standard values: init_bands(20,0).\n");
    init_bands(20,0);
  }
  
  if (bool_init_output == false) {
    print("init: Using standard values: init_output(180,90).\n");
    init_output(180, 90);
  }
  
  if (bool_init_rotations == false) {
    print("init: Using standard values: init_rotations_none() (disabled).\n");
    init_rotations_none();
  }
  
  if (bool_init_surface == false) {
    print("init: Using standard values: init_surface(10000).\n");
    init_surface(10000, true);
  }
  
  if (bool_init_translations == false) {
    print("init: Using standard values: init_translations(0) (disabled).\n");
    init_translations(0);
  }
  
  if (bool_init_slices == false) {
    print("init: Using standard values: init_slices(1) (disabled).\n");
    init_slices(1);
  }
  
  if (bool_init_noise == false) {
    print("init: Using standard values: init_noise(0) (disabled).\n");
    init_noise(0, ZERO);
  }
  
  if (bool_init_ocamcalib == false) {
    print("init: Using standard values: init_ocamcalib("") (disabled).\n");
  }   
  
  bool_init = true;
  
    
  // init some parameters and rotation stuff
  if (n_bands < 1) {
    n_bands = 1;
    print_warning("init", "n_bands must at least be 1");
  }
  if (n_bands_CG > n_bands) {
    print_warning("init", "n_bands_CG cannot be greater than n_bands");
    n_bands_CG = n_bands;
  } 
  
  // init stuff
  init_lib_parameters();
  
  print("init: Spherical Harmonics (surface): ");
  init_lib_surface();
  print("\n");
  
  print("init: Spherical Harmonics (output): ");
  init_lib_output();
  print("\n");

  print("init: Noise: ");
  init_lib_noise();
  print("\n");
  
  print("init: Rotations: ");
  for (int q=0; q<n_rot_par; q++) {
    init_lib_rotations(q);
  }
  print("\n");
  
  print("init: Translations: ");
  init_lib_translations();
  print("\n");
  
  print("init: Clebsch Gordan Matrices: ");
  init_lib_tm();
  init_lib_cg();
  print("\n");
  
  print("-------- Initialization Finished --------\n"); 
  
}

Shc::Shc() {
  
  bool_init_rotations              = false;
  bool_init_surface                = false;
  bool_init_output                 = false;
  bool_init_bands                  = false;
  bool_init_tolerances             = false;
  bool_init_files                  = false;
  bool_init_translations           = false;
  bool_init_slices                 = false;
  bool_init_noise                  = false;
  bool_init_ocamcalib              = false;
  bool_init                        = false;
  print_level                      = ALL;
  
  bool_feature_normalize           = false;
  feature_norm                     = NORM_2;
  bool_linearize_translations      = true;
  
  bool_measureOngoing              = false;
  bool_tick                        = false;
  
  rescale_factor                   = 2.0f;
  
  callback_filter_mat              = &filter_bilinear_mat;
  callback_filter_raw              = &filter_bilinear_raw;
  filter                           = BILINEAR;
  
  tangent_distance_mode            = OFF;
  tangent_distance_spring_constant = 0;
  
  noise_enable                     = false;
  noise_amplifier                  = 1.0f;
  
  bool_FFT                         = false;
  
  progress_current                 = -1;
  
  rand_generator.seed(time(NULL));  
  
}

