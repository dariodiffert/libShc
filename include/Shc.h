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

#ifndef AGTI_DD_SHC
#define AGTI_DD_SHC

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/FFT>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Shc_defs.h>

#include <kiss_fftr.h>

namespace shc {
  
  /**
   *  @brief Core class of the libShc.
   * 
   *  This class contains most functionality of the libShc.
   *  For an overview of all classes and their relations, see \ref classlayout.
   */  
  class Shc {

    friend class Shx;   
    
    private:
          
      // ------------------------------------------------------------------------------------------
      // ENUMS
      // ------------------------------------------------------------------------------------------
      
      enum e_surf {SURF_UNIFORM_SPHERE, SURF_UNIFORM_PANORAMA, SURF_CUSTOM};

      // ------------------------------------------------------------------------------------------
      // TYPEDEFS
      // ------------------------------------------------------------------------------------------
      
      typedef Eigen::Triplet<float> TD;
      typedef Eigen::Triplet<std::complex<float> > TC;
      typedef unsigned int uint;
      
      // ------------------------------------------------------------------------------------------
      // CLASSES
      // ------------------------------------------------------------------------------------------
      
      class Wfft {
  
        public:
          
          std::vector<unsigned char> buffer;
          int n;
          int inverse;
          
          Wfft();
          Wfft(int n, int inverse);
          void init(int n, int inverse);
          Wfft& operator=(const Wfft& src);
          Wfft(const Wfft& src);
          void compute(VectorReal& in, int offset, std::vector<kiss_fft_cpx>& out);
          void compute(std::vector<kiss_fft_cpx>& in, VectorReal& out, int offset);
          
      };
      
      // ------------------------------------------------------------------------------------------
      // STRUCTS
      // ------------------------------------------------------------------------------------------
      
      struct s_coor3 {
        float x;
        float y;
        float z;
        s_coor3() : x(0), y(0), z(0) {};
        s_coor3(float x, float y, float z) {this->x = x; this->y = y; this->z = z;};
        s_coor3(Coor3d in) {this->x = in(0); this->y = in(1); this->z = in(2);};
      };
      struct s_coor2 {
        float x;
        float y;
        s_coor2() : x(0), y(0) {};
        s_coor2(float x, float y) {this->x = x; this->y = y;};
        s_coor2(Coor2d in) {this->x = in(0); this->y = in(1);};
      };
      struct s_sphericalCoor {  
        float theta;
        float phi;
      };
      struct s_surf_unit {   
        s_sphericalCoor sphericalCoor; 
        s_coor3 coor;
        float weight;
      };
      struct s_surf {
        std::vector<s_surf_unit> unit;
        int n_angles;
        int n_sha;
        int l_max;
        VectorReal sh;
        VectorReal sh_weighted_full;
        VectorReal sh_weighted_hemi;
        e_surf type;
      };
      struct s_FFT {
        std::vector<Wfft> cfg_FFT;
        std::vector<Wfft> cfg_IFFT;
        std::vector<int> theta_index;
        std::vector<std::vector<kiss_fft_cpx>> store;
        VectorReal val_FFT;
        VectorReal val_IFFT;
      };
      struct s_IFFT_output {
        Wfft cfg;
        std::vector<kiss_fft_cpx> in;
        VectorReal val;
      };
      struct s_rotation_indices {
        int i;
        int j;
        int k;
      };
      struct s_compass_prepare_shpm {
        VecShpm cv;
        VecShpm ss;
        VecShpm rr_cur;
        VecShpm cv_cur;
        VecShpm ss_cur;
        int n_shpm;
        int n_mask;
        int n_total;
        std::vector< std::vector<int> > sh_max;
        std::vector< std::vector<int> > l_max;
        std::vector<bool> valid;
        std::vector<VecShpm> ss_td;
        std::vector<VecMatrixReal> L;
        std::vector<VecMatrixReal> LL;
      };
      struct s_translate_prepare {
        float w[2];
        float mu[2];
        std::vector<int> indices[2];
        std::vector<float> weights[2];
      };
      struct s_translate_unit {
        std::vector<std::vector<MatrixReal> > matrixDense;
        std::vector<std::vector<MatrixRealSparse> > matrixSparse;
        std::vector<std::vector<MatrixRealSparse> > matrixSparse_RM;
        std::vector<std::vector<MatrixRealSparse> > matrixSparse_RMN;
      };
      struct s_translate {
        std::vector<s_translate_unit> unit;
        VectorReal slices;
        VectorReal distances;
        int n_distances;
        int n_slices;
        int l_max;
        int n_points;
        e_contin contin;
        e_trans_type trans_type;
        float tolerance_translation;
      };
      struct s_surf_quickRef {
        int n_theta;
        float theta_stepSize;
        std::vector<float> theta;
        std::vector<int> theta_index;
        std::vector<int> psi_count;
      };
      // adapted from Scaramuzza (https://sites.google.com/site/scarabotix/ocamcalib-toolbox)
      struct s_ocam_model {
        double pol[1024];
        int length_pol;
        double invpol[1024];
        int length_invpol;
        double xc;
        double yc;
        double c;
        double d;
        double e;
        int width;
        int height;
        std::vector<s_coor2> mapped;
        bool invert_x_axis;
        bool invert_y_axis;
        bool invert_z_axis;
        MatrixReal mappingConversionY;
        MatrixReal mappingConversionX;
      };
      struct s_rot_Z {
        std::vector<VectorReal> v1;
        std::vector<VectorReal> v2;
      };
      struct s_rot_Y {
        VecMatrixReal M1;
        VecMatrixReal M2;
      };
      struct s_rot_X {
        VecMatrixReal M;
      };
      struct s_rot_X2 {
        VectorInt n[4];
        VecVectorInt p;
        VecVectorInt c[4];
        VecVectorInt r[4];
        VecMatrixReal M[4];
      };
      struct s_transform_unit {
        int l_max;
        e_matrix_type matrix_type;
        MatrixReal M;
        MatrixReal M_RM;
        MatrixReal M_RMN;
        MatrixRealSparse MS;
        MatrixRealSparse MS_RM;
        MatrixRealSparse MS_RMN;
      };
      struct s_transform {
        std::vector<s_transform_unit> unit;
      };
      struct s_rotation {
        VectorReal psi; 
        VectorReal theta;
        VectorReal phi;
        
        std::vector<s_rot_Z> rotate_Z;
        s_rot_Y rotate_YN;
        s_rot_Y rotate_YP;
        s_rot_X rotate_XN;
        s_rot_X rotate_XP;
        s_rot_X rotate_YNXN;
        
        std::vector<int> fast_indices_x;
        std::vector<int> fast_indices_y;
        std::vector<int> fast_indices_z;
        int n_fast_indices_x;
        int n_fast_indices_y;
        int n_fast_indices_z;
        int n_total_rotations;
        int l_max;
        int l_max_compass;
        float resolution;
        VectorInt td_active;
      }; 
      struct s_sh_index {
        int l;
        int m;
      };
      struct s_cgCoef {  
        MatrixReal cgIndices;
        VectorReal cgCoef;
        VectorReal column;
        VectorReal row;
        int n;
      };
      struct s_files {
        std::string surface;
        std::string translation;
        bool force_overwrite;
      };
      struct s_tolerances {
        float translation;
        float cg;
      };
 
      // ------------------------------------------------------------------------------------------
      // FUNCTIONS
      // ------------------------------------------------------------------------------------------
      
      void world2cam(s_coor2& coor2, s_coor3& coor3);
      void cam2world(s_coor2& coor2, s_coor3& coor3); 

      Surf coef2surf(int index);
      
      static MatrixReal intern_warp(MatrixReal& matrix, Coor3d translation, Xyz xyz, e_trans_type trans_type);
      
      static float intern_otsu(VectorReal& hist);
      static VectorReal intern_calcHistogram(MatrixReal& matrix, float theta_min, float theta_max, bool cumulative, bool panoramic);
      static float intern_histogramEqualization_calcEntry(VectorReal& hist, float val);
      MatrixReal intern_output2matrix(Surf& surf);
      
      void tangent_distance_calc_coef(Coef& result, Shpm& shpmE, Shpm& shpmP, MatrixReal& L, MatrixReal& LL);
      Shpm tangent_distance_calc(Shpm& shpmE, Shpm& shpmP);
      VecShpm tangent_distance_calc(VecShpm& shpmE, Shpm& shpmP);
      void tangent_distance_calc(Shpm& shpmE, Shpm& shpmP, Shpm& resultE, Shpm& resultP);
      void tangent_distance_prepare(Shpm& shpm, MatrixReal& L, MatrixReal& LL, VectorInt active);
      void tangent_distance_prepare(Shpm& shpm, MatrixReal& L, VectorInt active);
      Shpm intern_transform(Shpm& shpm, s_transform_unit& transform);
      
      bool intern_save_transform(s_transform& transform, std::string file);
      bool intern_load_transform(s_transform& transform, std::string file);
      
      int intern_add_transform(s_transform& transform, MatrixReal& matrix_transform, e_matrix_type matrix_type = DENSE, float tolerance = 0, int l_max = IGNORE);
      void intern_clear_transform(s_transform& transform);
            
      Shpm contin_convert(Shpm& shpm, e_contin c, int l_max = IGNORE);
      
      s_compass_prepare_shpm compass_prepare_shpm(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2);
      s_compass_prepare_shpm rotate_prepare_shpm(Shpm& shpm);
      s_compass_prepare_shpm rotate_prepare_shpm(VecShpm& vshpm);
      
      Surf intern_shpm2surf(Shpm& shpm);
      Surf intern_shpm2output(Shpm& shpm);
      
      void equalize_vector_bigger(Shpm& shpm1, Shpm shpm2, Coef& v1, Coef& v2);
      
      e_contin trifold(e_contin contin);
      
      template <typename X, typename Y>
      bool equal(X x, Y y) {return x == y;}

      template <typename X, typename Y, typename... Others>
      bool equal(X x, Y y, Others ... args) {return (x==y) && equal(y, args...);}

      bool check_init(bool is_init = false);
      bool check_init(Shpm& shpm);
      bool check_init(VecShpm& v_shpm);

      bool odd(int n);
      bool even(int n);
      bool contin_skip(e_contin contin, int l);

      void compass_calc_SS(VecShpm& result, s_compass_prepare_shpm& scp, int index_rot_par);
      void compass_calc_X(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par);
      void compass_calc_Y(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par);
      void compass_calc_Z(VecShpm& result, VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par);
      void compass_calc_CV(VecShpm& result, VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par);
      void rotate_calc_Z(VecShpm& in, int index, s_compass_prepare_shpm& scp, int index_rot_par);
      void rotate_calc_ZY(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par);
      void rotate_calc_YX(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par);
      void rotate_calc_XZ(VecShpm& in, s_compass_prepare_shpm& scp, int index_rot_par);
      
      VectorReal stepspace(float resolution, float max_val, bool negative);           
      void init_rotations_none();
      
      static void filter_bilinear_indices(float& y, float& x, float h, float w, int& y0, int& x0, int& y1, int& x1);
      static void filter_bilinear_getValue(MatrixReal& matrix, int y0, int x0, int y1, int x1, float& f00, float& f10, float& f01, float& f11);
      static void filter_bilinear_getValue(unsigned char* data, int width, int y0, int x0, int y1, int x1, float& f00, float& f10, float& f01, float& f11);
      static float filter_bilinear_calc(float f00, float f01, float f10, float f11, float y, float x);
      static void filter_nearest_indices(float y, float x, float h, float w, int& y0, int& x0);
      
      static float filter_bilinear_mat(MatrixReal& matrix, int h, int w, float y, float x);
      static float filter_bilinear_raw(unsigned char* data, int h, int w, float y, float x);
      static float filter_nearest_mat(MatrixReal& matrix, int h, int w, float y, float x);
      static float filter_nearest_raw(unsigned char* data, int h, int w, float y, float x);
      
      VecShpm intern_product(VecShpm& v_shpm1, VecShpm& v_shpm2);
      VecShpm intern_rotate(VecShpm& v_shpm, Xyz xyz);
      VecShpm intern_translate(VecShpm& v_shpm, float dist);
      VecShpm intern_translate_single(VecShpm& v_shpm, float dist);
      VecShpm intern_translate_dual(VecShpm& v_shpm, float dist);
      VecShpm intern_warp(VecShpm& v_shpm, Coor3d translation);
      
      std::vector<Coor3d> read_file_vec3(std::string file);
      
      static s_coor2 sphericalCoor2coor2(s_sphericalCoor sphericalCoor);
      static s_sphericalCoor coor22sphericalCoor(s_coor2 coor2);
      static s_coor3 sphericalCoor2coor3(s_sphericalCoor sphericalCoor);
      static s_sphericalCoor coor32sphericalCoor(s_coor3 coor3);
      
      std::vector<s_translate_prepare> calc_translation_prepare(s_coor3 trans, float slice, e_trans_type trans_type);
      std::vector<s_translate_prepare> calc_translation_prepare(float distance, float slice, e_trans_type trans_type);
      VecShpm calc_translation_column(std::vector<s_translate_prepare>& prepare, int l_max, int index, bool z_only, e_contin contin);
      VecMatrixReal calc_translation(float distance, float slice, e_trans_type trans_type, int l_max);
      MatrixReal calc_translation_simple(s_coor3 trans, e_trans_type trans_type, e_contin contin, int l_max);
      
      std::vector<int> get_nearest_indices(s_sphericalCoor& sphericalCoor);
      std::vector<float> get_nearest_weights(s_sphericalCoor& sphericalCoor, std::vector<int>& indices);
      float get_nearest_surf(Surf& surf, s_sphericalCoor& sphericalCoor);
      float get_nearest_surf(Surf& surf, std::vector<int>& indices, std::vector<float>& weight);
        
      s_surf_quickRef create_surface_quickRef(s_surf& surface);
      
      s_surf create_sphere(int n_points);
      s_surf create_sphere(int width, int height);
      
      int getClosestIndex(float angle, int index_max, VectorReal& angles, std::vector<int>& fastIndices);
      void xyz2rotationIndices(Xyz xyz, int index_rot_par, s_rotation_indices& indices);
      Xyz rotationIndices2xyz(int index_rot_par, s_rotation_indices& indices);
      
      float min3(float i1, float i2, float i3);
      float max3(float i1, float i2, float i3);
      void create_fast_indices(VectorReal angles, float resolution, s_rotation &rotation);
          
      MatrixReal reduce_matrix(MatrixReal& M, e_contin contin);
      
      bool calc_rotate_isEntryZero_X(int mi, int mj, int l_ninety);
      bool calc_rotate_isEntryZero_Y(int mi, int mj, int l_ninety);
      bool calc_rotate_isEntryZero_Z(int mi, int mj, int l_ninety);
      bool calc_rotate_isEntryZero(e_axis axis, int li, int lj, int mi, int mj, bool ninety = false);
      bool calc_rotate_isEntryZero(int li, int lj, int mi, int mj);
      bool calc_translate_isEntryZero(int mi, int mj);
      
      double calc_rotation_entry(MatrixRotation &R, MatrixReal &M, int l, int m, int n);
      double calc_rotation_entry_Z(double a,int mi, int mj);
      s_rot_Z create_rotate_Z(float rotation_angle, int l_max);
      s_rot_Y create_rotate_Y(MatrixReal& M, int l_max);
      s_rot_X create_rotate_X(MatrixReal& M, int l_max);
      s_rot_X2 create_rotate_X2(MatrixReal& M, int l_max);
      Coef compute_rotate_X2(s_rot_X2& rot, Coef& coef, int l_max);
      
      MatrixReal calc_rotate_simple(Xyz xyz, int l_max);
      MatrixReal calc_rotate_axis(float rotation_angle, e_axis axis, int l_max);
      double calc_rot_U(MatrixRotation &R, MatrixReal &M, int l, int m, int n);
      double calc_rot_V(MatrixRotation &R, MatrixReal &M, int l, int m, int n);
      double calc_rot_W(MatrixRotation &R, MatrixReal &M, int l, int m, int n);
      double calc_rot_P(MatrixRotation &R, MatrixReal &M, int i, int l, int a, int b);
      double calc_rot_u(int l, int m, int n);
      double calc_rot_v(int l, int m, int n);
      double calc_rot_w(int l, int m, int n);
        
      Xyz intern_compass(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2);
      float compass_evaluate_compare(s_compass_prepare_shpm& sc, int index_rot_par);
      s_rotation_indices compass_evaluate(int index_rot_par, s_compass_prepare_shpm& scp);
      
      float feature_difference(VectorReal& v1, VectorReal& v2);
      float calc_difference_ISE(Shpm& shpm1, Shpm& shpm2);
      float calc_difference_AS(Shpm& shpm1, Shpm& shpm2);
      float calc_difference_BS(Shpm& shpm1, Shpm& shpm2);
      float calc_difference_BS_dense(Shpm& shpm1, Shpm& shpm2);
      
      Shpm full2hemi(Shpm& shpm, e_contin contin);
      Shpm hemi2full(Shpm& shpm);
      
      Shpm intern_surf2shpm(Surf& surf, e_contin contin, int l_max, bool allow_noise, e_contin contin_result);
      Shpm intern_surf2shpm_hemi(Surf& surf, e_contin contin, int l_max);
      Shpm intern_surf2shpm_hemi_FFT(Surf& surf, e_contin contin, int l_max);
      Shpm intern_surf2shpm_full(Surf& surf, int l_max);
      Shpm intern_surf2shpm_full_FFT(Surf& surf, int l_max);
      Surf Shpm2surf(s_surf& surface, Shpm& shpm);
      
      void index2sh(int i, int& l, int& m);
      int sh2index(int l, int m);
      int sha2index(s_surf& surface, int l, int m, int a);
      int l2sh(int l, e_contin contin);

      bool init_lib_translations_load();
      bool init_lib_surface_load();
      
      void init_lib_angles(int n_points);
      void init_lib_output();
      void init_lib_rotations(int index_rot_par);
      void init_lib_tm();
      void init_lib_cg();
      void init_lib_parameters();
      void init_lib_translations();
      void init_lib_surface_FFT();
      void init_lib_surface_sh(s_surf& surface);
      void init_lib_surface();
      void init_lib_noise();
      
      double calc_sh_P(int l, int m, double x);
      double calc_sh_K(int l, int m);
      double calc_sh_Y(int l, int m, double theta, double phi);
         
      CoefBandwise shpm2bandwise(Shpm& shpm, int l_max = IGNORE);
      
      void cg_indexConversion(int l, int m, int l1, int m1, int l2, int m2, int& row, int& column);
      s_cgCoef cg_coefficients(int l, int m, int l1, int l2);
      MatrixRealSparse cg_createMatrix(int l1, int l2);
      MatrixRealSparse cg_createMatrixReal(int l1, int l2);
      MatrixRealSparse cg_createMatrixRealCoupling(int l1, int l2);
      MatrixComplexSparse createTransformationMatrix(int l);
      MatrixComplexSparse directSum(std::vector<MatrixComplexSparse>& in, int start, int stop);
      
      VectorReal calc_feature_ISE(Shpm& shpm);
      VectorReal calc_feature_amplitudespectrum(Shpm& shpm);
      VectorReal calc_feature_bispectrum(Shpm& shpm);
      VectorReal calc_feature_bispectrum_dense(Shpm& shpm);
      
      Shpm intern_product(Shpm& shpm1, Shpm& shpm2);
      
      bool save_lut_translate();
      bool load_lut_translate();
      bool save_lut_surf();    
      bool load_lut_surf();

      // ------------------------------------------------------------------------------------------
      // VARIABLES
      // ------------------------------------------------------------------------------------------
      
      s_ocam_model ocam_model;  

      s_FFT lut_FFT;
      s_IFFT_output lut_IFFT_output;
      int output_width;
      int output_height;
      
      bool bool_init_files;
      bool bool_init_tolerances;
      bool bool_init_rotations;
      bool bool_init_surface;
      bool bool_init_output;
      bool bool_init_bands;
      bool bool_init_noise;
      bool bool_init_translations;
      bool bool_init_slices;
      bool bool_init_ocamcalib;
      bool bool_init;
      
      bool bool_sphere_mode;
      bool bool_FFT;
      
      struct timespec measureStart;
      struct timespec measureTick;
      bool bool_measureOngoing;
      bool bool_tick;
      e_feature_norm feature_norm;
      bool bool_feature_normalize;
      bool bool_linearize_translations;
      e_print_level print_level;
      
      std::vector<MatrixComplexSparse> lut_tm;
      std::vector<std::vector<MatrixRealSparse> > lut_cg;
      std::vector<std::vector<MatrixRealSparse> > lut_cg_real;
      std::vector<std::vector<MatrixRealSparse> > lut_cg_realCoupling;
      
      std::vector<int> n_rotate_psi; 
      std::vector<int> n_rotate_theta;
      std::vector<int> n_rotate_phi;
      
      int n_rot_par;

      std::vector<s_rotation> v_rotation;

      int n_bands;
      int n_bands_CG;
      int n_coefficients;

      s_surf lut_surf; 
      
      s_tolerances tolerances;    
      s_files files;
      
      float rescale_factor;

      s_translate lut_translate;
      s_surf_quickRef lut_surf_quickRef;
      std::vector<s_sh_index> lut_index2sh;
      
      float progress_current;
      struct timespec progress_time;
      
      e_filter filter;
      float (*callback_filter_mat) (MatrixReal& matrix, int h, int w, float x, float y);
      float (*callback_filter_raw) (unsigned char* buffer, int h, int w, float x, float y);
      
      e_tangent_distance_mode tangent_distance_mode;
      float tangent_distance_spring_constant;
      s_transform tangent_distance;
      s_transform transforms;
      
      VectorReal noise_amplitudes;
      e_noise noise_type;
      int noise_samples;
      std::vector<Surf> noise_surf;
      Surf noise_mask;
      float noise_mask_sum;
      bool noise_enable;
      float noise_amplifier;
      
      std::default_random_engine rand_generator;
        
      
    protected:
      
      /**
       *  \brief Prints a message to the console
       *
       *  If the message is printed depends on the printing level, see Shc::set_print_level.
       *
       *  \param[in] message The message to print
       */
      void print(std::string message);    
      /** \overload
       *
       *  Prints a result of a computation to the console.
       *  For description, a prefix and suffix can be added.
       *
       *  \param[in] prefix Prefix text
       *  \param[in] result A floating number, commonly the result of some computation
       *  \param[in] suffix Suffix text
       */
      void print(std::string prefix, float result, std::string suffix);
      /** \overload
       *
       *  Calling this function, a progress number can be printed to the console via Shc::print_progress_start.
       *  The progress is set to 0%.
       */
      void print_progress_start();
      /** \overload
       *  
       *  Prints the current progress number to the console
       *
       *  \param[in] progress The progress in the range [0,1]
       */
      void print_progress_update(float progress);
      /** \overload
       *
       *  Stops printing a progress number to the console
       *
       *  Calling this function, the progress printing is stopped.
       *  As a consequence, Shc::print_progress_update cannot longer be called until a Shc::print_progress_start is called again.
       */
      void print_progress_stop();
      /** \overload
       *
       *  Prints a warning to the console
       *
       *  \param[in] function The function from which the warning originates
       *  \param[in] warning The warning message
       */
      void print_warning(std::string function, std::string warning);
      /** \overload
       *
       *  Prints a warning originating from a static function to the console
       *
       *  \param[in] function The static function from which the warning originates
       *  \param[in] warning The warning message
       */
      static void print_warning_static(std::string function, std::string warning);
            
    public:
    
      // initialization
      Shc();
      
      /**
       * \brief Initializes the number of bands (frequencies) used for computations in the basis of spherical harmonics
       * 
       * Translation and Clebsch-Gordan matrices (\ref translations, \ref products) can contain many values close to zero.
       * By setting the corresponding tolerance value, entries whose absolute values are small can also be set to zero.
       * Due to the increased sparsity of the matrices, this can reduce the computation times.
       * Note that by setting translation greater zero, translation matrices will be represented by sparse matrices instead of dense matrices.
       * This might actually increase the computation times if the tolerance value is too small.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] translation Threshold for translation matrices
       * \param[in] cg          Threshold for Clebsch-Gordan matrices
       */
      void init_tolerances(float translation, float cg);
      /**
       * \brief Initializes the file names in which precomputed values are saved
       * 
       * Surfaces and translation information (\ref surface, \ref translations) are precomputed.
       * To avoid the computation of these information at each initialization, these information can be stored in files.
       * After successfull initialization, the information are automatically stored in the specified files.
       * If the files exist at initialization, the information are loaded from the specified files.
       * If the loaded information differs from the used initialization settings, a warning is shown.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] file_surface Name of the file where the surface information is stored.
       * \param[in] file_translation Name of the file where the translation matrices are stored.
       * \param[in] force_overwrite If the loaded information differs from the used initialization settings, the initialization is redone with the current settings and afterwards saved to the specified files.
       *                            This setting avoids using accidentally old settings if these differ from the currently used settings.
       */
      void init_files(std::string file_surface, std::string file_translation, bool force_overwrite = true);
      /**
       * \brief Initializes the number of bands (frequencies) used for computations
       * 
       * Sets the maximal number of bands used for computations in the basis of spherical harmonics.
       * This setting affects <b>all</b> computations in the basis of spherical harmonics, e.g.
       * the Fourier transform, rotations, translations, etc.
       * Most functions have an additional parameter to use a smaller number of bands, however it is not possible to increase the number of bands after initialization.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] n_bands Maximal number of bands
       * \param[in] n_bands_CG Maximal number of bands used to compute Clebsch-Gordan matrices.
       *                       These matrices are necessary to compute products in the basis of spherical harmonics (\ref compass, \ref products).
       *                       Note that the calculation of the Clebsch-Gordan matrices is time consuming for higher bands.
       */
      void init_bands(int n_bands, int n_bands_CG = 0);
      /**
       * \brief Initializes the sets of rotation angles used to compute rotation matrices
       * 
       * This function can be used to create a set of rotation angles used to compute rotations in the basis of spherical harmonics (\ref rotations).
       * To implement rotations efficiently, rotations around all axes are computed via a set of Z-axis rotations and a fixed set of X-/Y-axis rotations around \f$ \pm 90^\circ \f$.
       * As a consequence, a shared set of \p rotation_angles is chosen to compute rotations around all axes.
       * Now for each axis a subset of these \p rotation_angles is chosen via \p indices_x, \p indices_y, and \p indices_z for each desired rotation.
       * Each subset is required to have at least one entry.
       * 
       * <b>Example (pseudo code):</b>
       * Assume that we want to allow rotations around the Z-axis in \f$ 90^\circ \f$ steps, i.e. \f$ rotation_angles = \{ 0^\circ, 90^\circ, 180^\circ, 270^\circ \} \f$.
       * Then we have to set \f$ indices_z = \{ 0,1,2,3 \} \f$.
       * Furthermore, we want to allow rotations around the Y-axis by either \f$ 0^\circ \f$ or \f$ 180^\circ \f$.
       * Then we have to set \f$ indices_y = \{ 0,2 \} \f$.
       * Finally, we are not interested in X-axis rotations, i.e. we only allow rotations by \f$ 0^\circ \f$.
       * Then we have to set \f$ indices_x = \{ 0 \} \f$.
       * 
       * <b>Coarse-to-fine:</b>
       * By calling this function multiple times, multiple sets of rotations are initialized.
       * Functions which compute rotations, e,g, Shc::rotate or Shc::compass, will use them consecutively.
       * This can be used to implement coarse-to-fine search strategies.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] rotation_angles A set of rotation angles
       * \param[in] indices_x A set of indices (corresponding to \p rotation_angles) for X-axis rotations
       * \param[in] indices_y A set of indices (corresponding to \p rotation_angles) for Y-axis rotations
       * \param[in] indices_z A set of indices (corresponding to \p rotation_angles) for Z-axis rotations
       * \param[in] l_max_compass The maximal number of bands used by Shc::compass. This parameter has to be smaller or equal to \p l_max. If set to zero, the compass skips this set of rotations
       * \param[in] l_max The maximal number of bands used
       */
      void init_rotations_custom(VectorReal rotation_angles, VectorInt indices_x, VectorInt indices_y, VectorInt indices_z, int l_max_compass = IGNORE, int l_max = IGNORE);
      /** \overload
       * 
       * This function can be used to create a predefined set of rotation angles.
       * This set is specified via the parameter \p resolution and \p range which are depicted in the following sketch:
       * 
       * \image html init_rotations_sphere.png
       * 
       * The term <i>sphere</i> refers to the idea that Z-axis rotations cover all directions.
       * In contrast, the X- and Y-axis use the same set of rotation parameters.
       */
      void init_rotations_sphere(float resolution, float range, int l_max_compass = IGNORE, int l_max = IGNORE);
      /** \overload
       * 
       * This function can be used to create a predefined set of rotation angles.
       * This set is specified via the parameter \p resolution and \p range which are depicted in the following sketch:
       * 
       * \image html init_rotations_cone.png
       * 
       * The term <i>cone</i> refers to the idea that rotations around all axes are limited to the same range.
       */
      void init_rotations_cone(float resolution, float range, int l_max_compass = IGNORE, int l_max = IGNORE);
      /**
       * \brief Initializes the surface used for the Fourier transform (equally distributed)
       * 
       * To compute the (discrete) Fourier transform, a set of sampling points has to be specified (\ref surface).
       * This function creates a set of sampling points equally distributed on the sphere.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] n_points The number of sampling points. Note that the actual number of sampling points can slightly differ.
       *                     The real number of sampling points can be obtained via Shc::get_surface_size.
       * \param[in] FFT If set to true, the <i>fast</i> Fourier transform is used.
       */
      void init_surface(int n_points, bool FFT = true);
      /**
       * \brief Initializes the surface used for the Fourier transform (grid)
       * 
       * To compute the (discrete) Fourier transform, a set of sampling points has to be specified (\ref surface).
       * This function creates a grid of \p width x \p height sampling points for the spherical coordinates.
       * Since the sampling points are not equally distributed, each sampling point has an (automatically computed) weighting value.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] width The number of sampling points in azimuth direction
       * \param[in] height The number of sampling points in altitude direction
       * \param[in] FFT If set to true, the <i>fast</i> Fourier transform is used.
       */
      void init_surface(int width, int height, bool FFT = true);
      /**
       * \brief Initializes the surface used for the Fourier transform (custom)
       * 
       * To compute the (discrete) Fourier transform, a set of sampling points has to be specified (\ref surface).
       * This function allows to define arbitrary sets of sampling points.
       * Each sampling point is required to have a weight.
       * This weight depends on the surface area of the sphere represented by each sampling point.
       * Each sampling point is represented by a triplet \p theta, \p phi, and \p weight.
       * Therefore all parameter vectors \p theta, \p phi, and \p weight have to be of same size.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] theta The \f$ \vartheta \f$ angles of the sampling points
       * \param[in] phi The \f$ \varphi \f$ angles of the sampling points
       * \param[in] weight The weights of the sampling points. The weights are automatically normalized.
       */
      void init_surface(VectorReal theta, VectorReal phi, VectorReal weight);
      /**
       * \brief Initializes the surface used for the Fourier transform (custom)
       * 
       * To compute the (discrete) Fourier transform, a set of sampling points has to be specified (\ref surface).
       * This function is similar to Shc::init_surface(VectorReal theta, VectorReal phi, VectorReal weight),
       * but reads the sampling points from a file (see \ref fileinput).
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] file_surface 
       */
      void init_surface(std::string file_surface);
      /**
       * \brief Initializes the dimension for images created via the inverse Fourier transform
       * 
       * For visualization or the approximation of \ref translations, the inverse Fourier transform has to be computed.
       * This function sets the dimension of the panoramic image created via the inverse Fourier transform.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] width Width of the resulting panoramic image
       * \param[in] height Height of the resulting panoramic image
       */
      void init_output(int width, int height);
      /**
       *  \brief Initializes a set of translation matrices
       * 
       *  Translations simulate the movement of a camera inside a unit sphere.
       *  Since \ref translations are expensive to compute in the basis of spherical harmonics, this function can be used to precompute a set of translations along the Z-axis.
       *  By combining rotations and translations, arbitrary translations can be computed.
       *  As a consequence, the rotations initialized via Shc::init_rotations should cover the complete space of 3D-rotations.
       *
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       *
       * \param[in] distances For these distances the movement of the camera inside a unit sphere is simulated 
       * \param[in] trans_type Specifies how the movement is simulated
       * \param[in] l_max The maximal number of bands used
       */
      void init_translations(VectorReal distances, e_trans_type trans_type = VISUAL, int l_max = IGNORE);
      /** \overload
       *
       * \param[in] n_distances To compute translations, the interval [0,1] is split into equal parts 
       * \param[in] trans_type Specifies how the movement is simulated
       * \param[in] l_max The maximal number of bands used
       */
      void init_translations(int n_distances, e_trans_type trans_type = VISUAL, int l_max = IGNORE);
      /** \brief Initializes multiple slices
       *
       * See \ref translations. 
       *
       * <b>Experimental:</b> This functionality has not been tested again since its implementation!
       * There is a high chance that bugs are present!
       *
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       *
       * \param[in] slices The radius of each slice
       */
      void init_slices(VectorReal slices);
      /** \overload
       *
       * \param[in] n_slices  To define a set of radii, the interval [0,1] is split into equal parts 
       */
      void init_slices(int n_slices);
      /**
       * \brief Initializes noise to fill in missing data in panoramic images
       * 
       * A set of \p n_samples panoramic images filled with \ref noise is created.
       * Allways when a panoramic image is Fourier transformed, one of these "noise images" is randomly chosen and used to fill in missing data.
       * The noise is created such that it has the given \p amplitude_spectrum.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] n_samples The number of noise samples created
       * \param[in] amplitude_spectrum The noise is created such that it has the given amplitude spectrum 
       */
      void init_noise(int n_samples, VectorReal amplitude_spectrum);
      /** \overload
       *
       * \param[in] n_samples The number of noise samples created
       * \param[in] noise Use a predefined type of noise 
       */
      void init_noise(int n_samples, e_noise noise);
      /**
       * \brief Initializes the parameters for unwrapping panoramic images
       * 
       * For unwrapping panoramic images, the OcamCalib can be used (\ref installation).
       * The calibration file, commonly called <i>calib_results.txt</i>, can be loaded to automatically unwrap images using Shx::load_shpm_ocamcalib, Shc::ocamcalib2matrix, Shc::ocamcalib2shpm, etc.
       * Make sure to use the correct \ref coordinatesystem such that the panoramic images can be unwrapped.
       * 
       * The parameters \p invert_x_axis, \p invert_y_axis, or \p invert_z_axis can be used to correct differences between the coordinate systems.
       * During our tests, it was necessary to invert all axes.
       * 
       * <b>Initialization:</b> This function is required to be called before the instance of the shc::Shc class is initialized (Shc::init).
       * If not used, a default set of parameters is used (see console output during initialization).
       * 
       * \param[in] file File name of the calibration results using the OcamCalib, commonly <i>calib_results.txt</i>.
       * \param[in] invert_x_axis Sets if the X-axis is inverted
       * \param[in] invert_y_axis Sets if the Y-axis is inverted
       * \param[in] invert_z_axis Sets if the Z-axis is inverted
       */
      void init_ocamcalib(std::string file, bool invert_x_axis = true, bool invert_y_axis = true, bool invert_z_axis = true);
      /**
       * \brief Finishes the initialization process
       * 
       * <b>Initialization:</b> Calling this function, the instance of the class shc::Shc is initialized using the previously defined settings (calls to functions of the form <i>init_*</i>).
       * After calling this function, no more calls to initialization functions are possible.
       */
      void init();

      /**
       * \brief Sets the print level
       * 
       * Information printed by the libShc, which originate from a specific shc::Shc instance, are either warnings or common messages.
       * By setting the print level, the printed information can be filtered.
       * 
       * \param[in] print_level Sets the print level
       */
      void set_print_level(e_print_level print_level);
       /**
       * \brief Gets the print level
       * 
       * See \ref Shc::set_print_level
       * 
       * \return Returns the current print level
       */
      e_print_level get_print_level();

      /**
       *  \brief Selects which parts of a panoramic image are filled with noise
       *
       *  If \ref noise has been initialized via Shc::init_noise, a mask can be defined to determine which pixel are filled with noise.
       *  The most intuitive way is to set the noise via a mask which represents a panoramic image.
       *  Each pixel can take values between [0,1], i.e. is represented by a grayscale panoramic image.  
       *  Pixel set to 1 are filled with noise, pixel set to 0 are filled with the original pixel information.
       *  Values inbetween linearly interpolate between noise and original pixel information.
       *  This can be useful to make more smooth transitions to areas which are filled with noise.
       *
       *  Internally, the mask is used to create a surface (shc::Surf) which represents the mask.
       *
       *  \param[in] mask The mask represents the areas in a panoramic image which 
       */
      void set_noise_mask(MatrixReal& mask);
      /** \overload
       *  
       *  Sets the noise mask as as if the panoramic image was captured with a camera which has an opening angle of \p angle in \f$ [0,2\pi] \f$.
       *  All pixel which are not visible by the camera are filled with noise.
       *
       *  \param[in] angle Opening angle of the camera objective in radians. 
       */
      void set_noise_mask(float angle);
      /** \overload
       *  
       *  Sets the noise mask directly for each sampling point on the \ref surface.
       *
       *  \param[in] surf The surface used as mask
       */
      void set_noise_mask(Surf& surf);
      /** \overload
       *
       *  \return The surface which defines for each sampling point if it is filled with \ref noise.
       */
      Surf get_noise_mask();
      /**
       *  \brief Enables or disables noise to fill in missing data in panoramic images
       *  
       *  If enabled, missing data in panoramic images is filled in with \ref noise.
       *  The type of noise can be set via the Shc::set_noise_mask functions.
       *  Note that noise has to be initialized via Shc::init_noise.
       *  In this case, noise is automatically enabled.
       *
       *  \param[in] enable Enables or disables the noise
       */
      void set_noise_enable(bool enable);
      /** \overload
       *
       *  \return Returns if \ref noise is enable or disables
       */
      bool get_noise_enable();
      
      /**
       *  \brief Multplies noise intensity by a factor
       *
       *  Note that \ref noise has to be initialized via Shc::init_noise.
       *  In this case, the amplifier value is set to 1.
       *  
       *  \param[in] amplifier Sets the noise amplifier (multiplier) to the given factor
       */
      void set_noise_amplifier(float amplifier);
      /** \overload
       *  
       *  \return Returns the noise amplifier (multiplier)
       */      
      float get_noise_amplifier();
      /** 
       *  \brief Enables/disables linearization for computing translations.
       *
       *  For precomputed \ref translations (Shc::init_translations, Shc::translate), only a discrete set of distances can be used.
       *  By enabling this option, translations are linearly interpolated.
       *  Note that this doubles the computation times.
       *
       *  \param[in] linearize Enables/disables linearization
       */ 
      void set_linearize_translations(bool linearize);
      /** \overload
       *  
       * \return Returns is linearization of translations is enabled/disabled.
       */
      bool get_linearize_translations();
      
      /**
       *  \brief Sets the filter used for interpolation
       *
       *  \param[in] filter Sets the filter used to interpolate between pixel (subpixel)
       */
      void set_filter(e_filter filter);
      /** \overload
       *  
       * \return Returns the filter used to interpolate between pixel
       */
      e_filter get_filter();
      
      /** 
       *  \brief Enables / disables normalization for feature comparisons
       *
       *  The difference between two features \f$ u,v \f$ (i.e. Fourier coefficient vectors, amplitude spectra, etc.) is computed either as
       *  \f$ \| u - v \| \f$ (normal) or \f$ \frac{\| u - v \|}{\|u\| + \|v\|} \f$ (normalized).
       *  The norm used can be set via Shc::set_feature_norm.
       *  This setting influences the functions Shc::compass and Shc::get_feature_difference.
       *
       *  \param[in] feature_normalize Enables/disables the normalization during feature comparisons
       */
      void set_feature_normalize(bool feature_normalize);
      /** \overload
       *
       *  \return Returns if normalization during feature comparisons is enabled or disabled
       */
      bool get_feature_normalize();

      /** 
       *  \brief Sets the norm used for feature comparisons
       *
       *  The difference between two features \f$ u,v \f$ (i.e. Fourier coefficient vectors, amplitude spectra, etc.) is computed either as
       *  \f$ \| u - v \| \f$ (normal) or \f$ \frac{\| u - v \|}{\|u\| + \|v\|} \f$ (normalized).
       *  This function can be used to change the norm.
       *  To enable/disable normalization use Shc::set_feature_normalize.
       *  This setting influences the functions Shc::compass and Shc::get_feature_difference.
       *
       *  \param[in] feature_norm Sets the norm used for feature comparisons
       */
      void set_feature_norm(e_feature_norm feature_norm);
      /** \overload
       *
       *  \return Returns the norm used for feature comparisons
       */
      e_feature_norm get_feature_norm();
      
      /** 
       *  @brief Prints a short description to the console
       *
       *  Prints a short description of the passed instance to the console.
       *  Note that this function is not affected by the setting of Shc::set_print_level.
       *
       *  @param[in] xyz Rotation to print
       */  
      static void print(Xyz xyz);
      /** \overload
       *  
       *  @param[in] r Rotation matrix to print
       */  
      static void print(MatrixRotation r);
      /** \overload
       *  
       *  @param[in] axr Rotation to print
       */  
      static void print(Axr axr);
      /** \overload
       *  
       *  @param[in] shpm Prints a description of the Fourier coefficient vector
       */  
      static void print(Shpm shpm);
      /** \overload
       *  
       *  @param[in] pointcloud Prints a description of the pointcloud
       */  
      static void print(Pointcloud pointcloud);
       
      /** 
       *  @brief Converts angles from radians to degrees
       *  @param[in] rad angle in radians
       *  @return Returns angle in degrees
       */  
      static inline float r2d(float rad) {return rad*180.0/M_PI;}
      /** 
       *  @brief Converts angles from degrees to radians
       *  @param[in] deg Angle in degrees
       *  @return Returns angle in radians
       */  
      static inline float d2r(float deg) {return deg*M_PI/180.0;}
          
      /** 
       *  @brief Starts a time measurement
       *  @return Returns true if a new measurement could be started. If false, Shc::tock() needs to be called first to stop the ongoing time measurement.
       */  
      bool tick();
      /** 
       *  @brief Stops a time measurements
       *  @return Returns the time in milliseconds passed since calling Shc::tick(). If -1 is returned, no time measurement was started before.
       */  
      float tock();
      
      /** 
       *  @brief Computes the rotation angle of the given rotation
       *
       *  @param[in] xyz Rotation of which the rotation angle is computed
       *  @return The rotation angle in radians
       */  
      static float rotation_angle(Xyz xyz);
      /** 
       *  @brief Normalizes the rotation angle
       *
       *  @param[in] angle Rotation angle
       *  @return The rotation angle \p angle normalized to the interval \f$ [0, 2\pi) \f$
       */ 
      static float angular_normalization(float angle);
      
      /** 
       *  @brief Converts rotation: XYZ -> R
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] xyz Rotation
       *  @return Converted rotation
       */ 
      static MatrixRotation xyz2r(Xyz xyz);
      /** 
       *  @brief Converts rotation: R -> XYZ
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] r Rotation
       *  @return Converted rotation
       */ 
      static Xyz r2xyz(MatrixRotation r);
      /** 
       *  @brief Converts rotation: AXR -> R
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] axr Rotation
       *  @return Converted rotation
       */ 
      static MatrixRotation axr2r(Axr axr);
      /** 
       *  @brief Converts rotation: R -> AXR
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] r Rotation
       *  @return Converted rotation
       */ 
      static Axr r2axr(MatrixRotation r);
      /** 
       *  @brief Converts rotation: AXR -> XYZ
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] axr Rotation
       *  @return Converted rotation
       */ 
      static Xyz axr2xyz(Axr axr);
      /** 
       *  @brief Converts rotation: XYZ -> AXR
       *
       *  For details between the different representations for rotations, see \ref rotationparameters.
       *
       *  @param[in] xyz Rotation
       *  @return Converted rotation
       */ 
      static Axr xyz2axr(Xyz xyz);
      
      /** 
       *  @brief Represents a vector via its spherical coordinate
       *
       *  Represents a normalized vector by its spherical coordinates \f$ (\vartheta, \varphi) \f$ (physics convention).
       *  If the vector is a null vector, the spherical coordinate \f$ \pi/2, 0 \f$ is returned.
       *
       *  @param[in] vector An arbitrary vector, it is automatically normalized
       *  @param[out] theta The \f$ \vartheta \f$ angle which represents \p vector in spherical coordinates
       *  @param[out] phi The \f$ \varphi \f$ angle which represents \p vector in spherical coordinates
       */ 
      static void vector2sphericalCoordinate(Coor3d vector, float& theta, float& phi);
      /** 
       *  @brief Represents a spherical coordinate via a vector
       *
       *  Represents a spherical coordinate (physics convention) via a normalized vector
       *
       *  @param[in] theta The angle \f$ \vartheta \f$
       *  @param[in] phi The angle \f$ \varphi \f$
       *  @return A normalized vector pointing in direction of the spherical coordinate
       */ 
      static Coor3d sphericalCoordinate2vector(float theta, float phi);

      /**
       *  \brief Creates a tilt matrix 
       *
       *  Creates a rotation matrix which "tilts" the image (rotation around X/Y-axes).
       *  Using the example from \ref coordinatesystem, a tilt of the camera would appear if a wheel of the robot is elevated.  
       *
       *  \param[in] tilt_dir Direction in which the tilt occurs (the rotation axis is orthogonal to this angle)
       *  \param[in] tilt_angle Rotation angle
       *  \return The tilt matrix
       */      
      static MatrixRotation tilt2r(float tilt_dir, float tilt_angle);
      /** \overload
       *  
       *  See \ref rotationparameters.
       */
      static Xyz tilt2xyz(float tilt_dir, float tilt_angle);
      /** \overload
       *  
       *  See \ref rotationparameters.
       */
       static Axr tilt2axr(float tilt_dir, float tilt_angle);
      
      /** 
       *  @brief Transposes (inverts) a rotation
       *
       *  Computes the inverse of a rotation; for a rotation matrix this is the transposed.
       *  This function is overloaded to work with the various \ref rotationparameters.
       *
       *  @param[in] r The rotation
       *  @return The inverted rotation
       */ 
      static MatrixRotation transpose(MatrixRotation r);
      /** \overload
       */
      static Xyz transpose(Xyz xyz);
      /** \overload
       */
      static Axr transpose(Axr axr);
      
      
      /** 
       *  @brief Computes the product (concatenation) of two rotations.
       *
       *  This function is overloaded to work with the various \ref rotationparameters.
       *  @param[in] r1 The first rotation matrix
       *  @param[in] r2 The second rotation matrix
       *  @return Returns the product \p r = \p r1 * \p r2
      */ 
      static MatrixRotation product(MatrixRotation r1, MatrixRotation r2);
      /** \overload
       */
      static Xyz product(Xyz xyz1, Xyz xyz2);
      /** \overload
       */
      static Axr product(Axr axr1, Axr axr2);
      
      /**
       *  \brief Generates a linearly spaced vector
       *  
       *  Creates a linearly spaced vector on the interval [a,b] with \p steps points.
       *  The space between all points is similar.
       *
       *  \param[in] a Lower limited
       *  \param[in] b Upper limited
       *  \param[in] steps Number of points
       *  \return The linearly spaced vector as specified above
       */
      static VectorReal linspace(float a, float b, int steps);
      /**
       *  \brief Computes a trapezoidal numerical integration
       *  
       *  For a function \f$ f(x) = y \f$ with descrete sampling points,
       *  the integral is computed via trapezoidal numerical integration.
       *  The input vectors \p x and \p y are required to have the same size.
       *
       *  \param[in] x Sampling points of the function \f$ f \f$
       *  \param[in] y Values \f$ y = f(x) \f$ at the given sampling points
       *  \return The approximated integral of \f$ f \f$
       */
      static float trapz(VectorReal& x, VectorReal& y);

      /** 
       *  \brief Computes the area on the sphere which is visible in two rotationally misaligned hemispherical images
       *
       *
       *  \param[in] xyz1 Orientation of the first hemispherical panoramic image
       *  \param[in] xyz2 Orientation of the second hemispherical panoramic image
       *  \return The percentage of the area on the sphere which is visible in both hemispherical panoramic images
       */
      static float hemispherical_intersection(Xyz xyz1, Xyz xyz2);
      /** \overload
       *  
       *  \param[in] xyz Orientations for a set of hemispherical panoramic images
       *  \return The percentage of the area on the sphere which is visible in all hemispherical panoramic images
       */
      static float hemispherical_intersection(std::vector<Xyz> xyz);
      
      /**
       *  \brief Computes the smallest rotation angle necessary to transfer two rotations into each other
       *
       *  Given two rotations \f$ R_1, R_2 \f$, this function computes the rotation angle of
       *  \f$ R_1^T R_2 \f$, see shc::Axr.
       *  The angle is in the interval \f$ [0,\pi] \f$.
       *  This function is overloaded to work with the various \ref rotationparameters.
       *
       *  \param[in] r1 First rotation matrix
       *  \param[in] r2 Second rotation matrix
       *  \return The rotation angle between the rotations \p r1 and \p r2.
       */
      static float angular_difference(MatrixRotation r1, MatrixRotation r2);
      /// /overload
      static float angular_difference(Xyz xyz1, Xyz xyz2);
      /// /overload
      static float angular_difference(Axr axr1, Axr axr2);
      /** /overload
       *  
       *  Note: This function is limited to 2D rotations.
       */
      static float angular_difference(float angle1, float angle2);
      
      /** 
       *  \brief Returns the surface information
       *  
       *  After initialization via Shc::init_surface, this function can be used to obtain the \ref surface information
       *  The surface information is returned as a triple \f$ (\vartheta, \varphi, w) \f$ representing the location in spherical coordinates and the weight.
       *
       *  \param[out] theta The altitude angles \f$ \vartheta \f$ (spherical coordinates) of the sampling points
       *  \param[out] phi The azimuth angles \f$ \vartheta \f$ (spherical coordinates) of the sampling points
       *  \param[out] weight The weights of the sampling points
       */
      void get_surface_data(VectorReal& theta, VectorReal& phi, VectorReal& weight);
      /** 
       *  \brief Returns the number of sampling points
       *  
       *  After initialization via Shc::init_surface, this function can be used to obtain the number of sampling points used.
       *
       *  \return The number of sampling points
       */
      int get_surface_size();
      /** 
       *  \brief For a given spherical coordinate, this function returns the closest sampling point
       *  
       *  After initialization via Shc::init_surface, this function can be used to find the sampling point closest to a given spherical coordinate \f$ (\vartheta, \varphi) \f$.
       *
       *  \param[in] theta Altitude angle \f$ \vartheta \f$ (spherical coordinates)
       *  \param[in] phi Azimuth angle \f$ \varphi \f$ (spherical coordinates)
       *  \return Index of the closest sampling point, if required, use Shc::get_surface_data to get more detailed information
       */      
      int get_surface_nearest(float theta, float phi);
      /** \overload
       *
       *  \param[in] coor Coordinate (automatically normalized) of a point on the sphere
       *  \return Index of the closest sampling point, if required, use Shc::get_surface_data to get more detailed information
       */      
      int get_surface_nearest(Coor3d coor);
      /**
       *  \brief Normalizes the surface
       *
       *  \param[in] surf The shc::Surface to normalize
       *  \return A surface where all entries are normalized to he interval [0,1]
       */
      Surf normalize(Surf& surf);
      
      /**
       *  \brief Conversion: matrix > surf
       *
       *  \param[in] matrix Matrix (i.e. a panoramic image) to convert
       *  \return The resulting surface
       */
      Surf matrix2surf(MatrixReal& matrix);
      /** \overload
       *
       *  For higher performance, this function can be used to directly convert a matrix into a surface.
       *
       *  \param[in] data Pointer to a buffer holding the panoramic image information
       *  \param[in] width Width of the panoramic image
       *  \param[in] height Height of the panoramic image
       *  \return The resulting surface
       */
      Surf matrix2surf(unsigned char* data, int width, int height);
      /**
       *  \brief Conversion: shpm > matrix
       *
       *  \param[in] shpm shc::Shpm to convert
       *  \return The resulting matrix (inverse Fourier transform)
       */      
      MatrixReal shpm2matrix(Shpm& shpm);
      /**
       *  \brief Conversion: matrix > shpm
       *
       *  \param[in] matrix The panoramic image to Fourier transform
       *  \param[in] contin \ref hemisphericalcontinuations used for the Fourier transform
       *  \param[in] l_max Maximal number of bands used for the Fourier transform
       *  \return The resulting shc::Shpm (Fourier transform)
       */      
      Shpm matrix2shpm(MatrixReal& matrix, e_contin contin = FULL, int l_max = IGNORE);
      /** \overload
       *
       *  For higher performance, this function can be used to directly operate on a buffer.
       *
       *  \param[in] data Pointer to a buffer holding the panoramic image information
       *  \param[in] width Width of the panoramic image
       *  \param[in] height Height of the panoramic image
       *  \param[in] contin \ref hemisphericalcontinuations used for the Fourier transform
       *  \param[in] l_max Maximal number of bands used for the Fourier transform
       *  \return The resulting shc::Shpm (Fourier transform)
       */      
      Shpm matrix2shpm(unsigned char* data, int width, int height, e_contin contin = FULL, int l_max = IGNORE);
      
      /** 
       *  \brief Resizes the matrix to the given dimensions
       *  
       *  Each pixel in the resized matrix is computed using bilinear filtering.
       *  For high quality resizing of images, use appropriate software.
       *   
       *  \param[in] matrix Input matrix to resize
       *  \param[in] width Width of the resulting matrix
       *  \param[in] height height of the resulting matrix
       *  \return Returns the resized matrix with the given dimensions.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal resize(MatrixReal& matrix, int width, int height);
      /** 
       *  \brief Computes the pixel-wise product between two matrices.
       *  
       *  The pixel-wise product is commutative, i.e. the order of the input arguments is not important.
       *  Note that the matrices \p matrix1 and \p matrix2 are required to have the same dimensions.
       * 
       *  \param[in] matrix1 First input matrix
       *  \param[in] matrix2 Second input matrix
       *  \return Returns the result as a matrix with the same dimensions as the input matrices.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal product(MatrixReal& matrix1, MatrixReal& matrix2);
      /** 
       *  \brief Computes the rotation of a panoramic image (represented by a matrix)
       *  
       *  Applies a rotation to a panoramic image, represented by a matrix, according to the \ref coordinatesystem.
       *  For the spatial domain, this is equivalent to \ref rotations of spherical harmonics in the frequency domain.
       * 
       *  \param[in] matrix Input matrix which represents a panoramic image
       *  \param[in] xyz Specifies the rotation, see shc::Xyz.
       *  \return Returns the rotated panoramic image.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal rotate(MatrixReal& matrix, Xyz xyz);
      /** 
       *  \brief Warps/translates a panoramic image (represented by a matrix)
       *  
       *  Applies a translation (also called warping) to a panoramic image, represented by a matrix, according to the \ref coordinatesystem.
       *  For the spatial domain, this is equivalent to \ref translations of spherical harmonics in the frequency domain.
       *  
       *  \param[in] matrix Input matrix which represents a panoramic image
       *  \param[in] translation Target coordinate
       *  \param[in] trans_type Defines how the translation is interpreted (see \ref translations) 
       *  \return Returns the warped panoramic image.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal warp(MatrixReal& matrix, Coor3d translation, e_trans_type trans_type = VISUAL);
      /** 
       *  \brief Normalizes the matrix
       *  
       *  All entries in the matrix are normalized to the interval [0,1], 
       *  
       *  \param[in] matrix Input matrix
       *  \return Returns the normalized matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal normalize(MatrixReal& matrix);
      /** 
       *  \brief Convolutes the matrix with a given kernel
       *  
       *  This function convolutes a matrix with another given kernel matrix,
       *  e.g. to compute edge filtered images.
       *  At borders, the image is wrapped.
       *  
       *  \param[in] matrix Input matrix
       *  \param[in] kernel Kernel matrix. The kernel matrix must have odd dimensions.
       *  \return Returns the convoluted matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal convolution(MatrixReal& matrix, MatrixReal &kernel);
      /** \overload 
       * 
       *  This functions is optimized for 3x3 kernel matrices
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] kernel A 3x3 kernel matrix
       *  \return Returns the convoluted matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       * 
       */
      static MatrixReal convolution(MatrixReal& matrix, Eigen::Matrix3f& kernel);
      /** \overload 
       * 
       *  This functions provides a set of predefined convolutions.
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] convolution_type Specifies which kernel matrix should be used for the convolution
       *  \return Returns the convoluted matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       * 
       */
      static MatrixReal convolution(MatrixReal& matrix, e_convolution_type convolution_type);
      /**
       * \brief Applies the Sobel filter
       * 
       * The matrix is convoluted with the Sobel Filter, compare convolution(MatrixReal& matrix, e_convolution_type convolution_type)
       * 
       *  \param[in] matrix Input matrix
       *  \return Returns the Sobel filtered matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal sobel(MatrixReal& matrix);
      /**
       * \brief Applies histogram equalization
       * 
       *  The matrix is histogram equalized. 
       *  By setting the parameters \p bins_x and \p bins_y to values greater one, the image is divided into that many bins to perform (bilinearly interpolated) local histogram equalization.
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] bins_x Number of bins in x-direction
       *  \param[in] bins_y Number of bins in y-direction
       *  \param[in] panoramic The matrix represents a panoramic image, i.e. distortions towards the poles are factored in during the computation of the requreid histogram
       *  \return Returns the histogram equalized matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       */
      static MatrixReal histogram_equalization(MatrixReal& matrix, int bins_x = 1, int bins_y = 1, bool panoramic = true);
      /**
       *  \brief Applies thresholding to the matrix
       *
       *  All entries in the matrix which are smaller or greater than the given \p threshold are set to 0 or 1, respectively.
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] threshold The threshold value
       *  \return Returns the thresholded matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       *
       */
      static MatrixReal thresholding(MatrixReal& matrix, float threshold);
      /**
       *  \brief Applies automatic thresholding to the matrix
       *
       *  All entries in the matrix which are smaller or greater than the some threshold are set to 0 or 1, respectively.
       *  The threshold is computed via the Otsu algorithm or a variant of it, for more details see <i>Differt and Möller (2015)</i> in \ref literature.
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] lambda If set to zero, the Otsu algorithm is used to determine an optimal threshold.
       *                    For values greater zero, the threshold is readjusted as described in <i>Differt and Möller (2015)</i>.
       *  \param[in] panoramic The matrix represents a panoramic image, i.e. distortions towards the poles are factored in during the computation of the requreid histogram
       *  \return Returns the thresholded matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       *
       */
      static MatrixReal thresholding_otsu(MatrixReal& matrix, float lambda = 0, bool panoramic = true);
      /**
       *  \brief Applies directional blur to the panoramic image
       *
       *  The matrix (which represents a panoramic image) is blurred as if the camera would be rotated around the given axis.
       *  The blur is computed using a Gaussian with the parameter \p sigma, i.e. \f$ \frac{1}{\sqrt{2 \pi s^2}}  e^{-\frac{(x-c)^2}{2 s^2}} \f$
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] axis Specifies the axis around which the camera is rotated to create blur
       *  \param[in] sigma Specifies the standard deviation of the Gaussian used to blur the image
       *  \param[in] inverse Experimental (and naive) support for deblurring, better avoid it
       *  \return Returns the blurred matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       *
       */
      static MatrixReal blur(MatrixReal& matrix, Coor3d axis, float sigma, bool inverse);     
      
      /**
       *  \brief Projects a panoramic image onto a plane
       *
       *  The matrix (which represents a panoramic image) is projected on a plane perpendicular to the Z-axis.
       *  Note that this function is not optimized is mainly used for debbuging purposes. 
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] angle The maximal viewing angle projected onto the plane, limited to \f$ (0^\circ,180^\circ) \f$.
       *  \param[in] size The size in pixel of the resulting square matrix
       *  \return Returns the projected matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       *
       */
      static MatrixReal project_plane(MatrixReal& matrix, float angle, int size);
      /**
       *  \brief Creates an artificial fish-eye image from a panoramic image
       *
       *  The matrix (which represents a panoramic image) is projected into a typical fish-eye image.
       *  The viewing angle increases constantly with the radius of the fish-eye image.
       * 
       *  \param[in] matrix Input matrix
       *  \param[in] angle The maximal viewing angle of the fish-eye image.
       *                   The angle is limited to \f$ (0^\circ, 360^\circ) \f$, common fisheye lenses have have an opening angle of \f$ 120^\circ \f$ to \f$ 220^\circ\f$.
       *  \param[in] size The size in pixel of the resulting square matrix
       *  \return Returns the projected matrix.
       *          If the operation was not successfull, an empty matrix is returned.
       *
       */
      static MatrixReal project_fisheye(MatrixReal& matrix, float angle, int size);
      
      /**
       *  \brief Loads a set of points (3D coordinates) from a file
       *
       *  See \ref fileinput for the required text file format.
       *
       *  \param[in] file_pointcloud Filename of the text file in which the points are stored
       *  \return The set of all points
       */
      Pointcloud load_pointcloud(std::string file_pointcloud);
      /** \overload
       *
       *  Loads a set of points (3D coordinates) from a file, centers them (mean = 0), and scales them such that the distance of the most distant point is 1.
       */
      Pointcloud load_pointcloud_auto(std::string file_pointcloud);
      /** 
       *  \brief Rotates a set of points (3D coordinates)
       *
       *  \param[in] pointcloud A set of points
       *  \param[in] rotation Rotation applied to all points
       *  \return The rotated points
       */
      Pointcloud rotate(Pointcloud& pointcloud, Xyz rotation);
      /** 
       *  \brief Translates (warps) a set of points (3D coordinates)
       *
       *  \param[in] pointcloud A set of points
       *  \param[in] translation The translation applied to all points
       *  \param[in] scale Scales all points (relative to the origin) by the given factor
       *  \return The translated and scaled points
       */
      Pointcloud warp(Pointcloud& pointcloud, Coor3d translation, float scale = 1);
      
      /** 
       *  \brief \brief Conversion: pointcloud > shpm
       *
       *  Maps a set of points onto the surface of a unit sphere (histogram) and computes the Fourier transform vector
       *  If multiple slices are used (Shc::init_slices), the points are (linearly weighted) mapped to the two closest slices.
       *
       *  \param[in] pointcloud A set of points
       *  \param[in] contin The \ref hemisphericalcontinuations used
       *  \param[in] l_max The maximal number of bands used
       *  \return A vector of shc::Shpm instances, one for each slice
       */
      VecShpm pointcloud2shpm(Pointcloud& pointcloud, e_contin contin = FULL, int l_max = IGNORE);
      /**
       *  \brief Creates a rotation matrix in the basis of spherical harmonics
       *
       *  Creates an explicit rotation matrix (real Wigner-D matrix) which can be used to rotate a function in the basis of real spherical harmonics (see \ref rotations).
       *  The resulting square matrix is of dimension \f$ l^2 \times l^2 \f$, where \f$ l \f$ is the number of bands used.
       *
       *  \param[in] xyz The rotation, see \ref rotationparameters
       * \param[in] l_max The maximal number of bands used
       *  \return The resulting rotation matrix
       */
      MatrixReal create_matrix_rotation(Xyz xyz, int l_max = IGNORE);
      /**
       *  \brief Creates a translation matrix in the basis of spherical harmonics
       *
       *  Creates an explicit translation matrix which can be used to translate a function in the basis of real spherical harmonics (see \ref translations).
       *  The resulting square matrix is of dimension \f$ l^2 \times l^2 \f$, where \f$ l \f$ is the number of bands used.
       *
       *  \param[in] translation The goal location
       *  \param[in] trans_type Specifies how the movement is simulated
       *  \param[in] contin The \ref hemisphericalcontinuations used
       *  \param[in] l_max The maximal number of bands used
       *  \return The resulting translation matrix
       */
      MatrixReal create_matrix_warp(Coor3d translation, e_trans_type trans_type = VISUAL, e_contin contin = FULL, int l_max = IGNORE);
      
      /** 
       *  \brief Enables / disables the tangent distance and sets its mode and spring constant
       *
       *  The tangent distance is by default disabled.
       *  As stated in the original paper (\ref literature), the \p spring_constant can be set to decrease numerical problems throughout the computation.
       *
       *  \param[in] mode Enables / disables the tangent distance and sets the mode (one-sided, two-sided)
       *  \param[in] spring_constant Value of the spring constant
       */
      void set_tangent_distance(e_tangent_distance_mode mode, float spring_constant = 0);
      /** \overload
       *  
       *  \param[out] mode Current mode set
       *  \param[out] spring_constant Current value of the spring constant
       */
      void get_tangent_distance(e_tangent_distance_mode& mode, float& spring_constant);
      
      /**
       *  \brief Adds an arbitrary tangent distance transformation matrix
       *
       *  This function can be used to add transformation matrices used to compute tangent distance transformation matrices, compare Shc::add_transform. 
       */
      int add_tangent_distance(MatrixReal& transform, e_matrix_type matrix_type = DENSE, float tolerance = 0, int l_max = IGNORE);
      /**
       *  \brief Creates a rotation matrix which is then added as tangent distance transformation matrix
       *  
       *  The rotation matrix is computed as \f$ \frac{1}{2} \left(\mathbf{R}_{v,\alpha} + \mathbf{R}_{v,-\alpha}\right) \f$,
       *  where \f$ v \f$ is the rotation axis and \f$ \alpha \f$ the rotation angle.  
       *
       *  \param[in] axis The rotation axis
       *  \param[in] angle The rotation angle
       *  \param[in] l_max The maximal number of bands used
       *  \return If successfull, the index of the added tangent distance transformation is returned. Otherwise, -1 is returned.
       */
      int add_tangent_distance_rotation(e_axis axis, float angle, int l_max = IGNORE);
      /**
       *  \brief Creates a translation matrix which is then added as tangent distance transformation matrix
       *  
       *  The translation matrix is computed as \f$ \frac{1}{2} \left(\mathbf{R}_{v,d} + \mathbf{R}_{v,-d}\right) \f$,
       *  where \f$ v \f$ is the axis of translation (direction) and \f$ d \f$ the distance.  
       *
       *  \param[in] axis The translation axis
       *  \param[in] dist The distance
       *  \param[in] contin \ref hemisphericalcontinuations used
       *  \param[in] l_max The maximal number of bands used
       *  \return If successfull, the index of the added tangent distance transformation is returned. Otherwise, -1 is returned.
       */
      int add_tangent_distance_translation(e_axis axis, float dist, e_contin contin = FULL, int l_max = IGNORE);
      /**
       *  \brief Returns the number of tangent distance transformation matrices loaded
       *
       *  \return The tangent distance transformation counter (transformations successfully added via Shc::add_tangent_distance and Shc::load_tangent_distance)
       */
      int get_tangent_distance_size();
      /**
       *  \brief Clears all tangent distance transformations, the tangent distance transformation counter is set to zero.
       */
      void clear_tangent_distance();
      /**
       *  \brief Configures for each stage of the visual compass which tangent distance transformations are applied
       *
       *  Due to the computational costs of using the tangent distance during the visual compass,  
       *  it can be specified for each stage of rotations (initialized via Shc::init_rotations) which tangent distance transformations are used.
       *  
       *  \param[in] index The index of the rotation stage (indices start at 0)
       *  \param[in] td_active For each entry of the vector, the tangent distance transformation with the same index is enabled.
       *                       All other tangent distance transformations are disabled for that rotation stage.
       */
      void compass_configure_tangent_distance(int index, VectorInt td_active);
      
      /**
       *  \brief Adds an arbitrary transformation matrix
       *
       *  This function can be used to add an arbitrary transformation matrix of size \f$ l^2 \times l^2 \f$, where \f$ l \f$ is the maximal number of bands initialized via Shc::init_bands.
       *  The matrix is automatically preprocessed such that it can be used for all \ref hemisphericalcontinuations.
       *  Transformation matrices can for example be created using Shc::create_matrix_rotation and Shc::create_matrix_warp.
       *  To avoid creating transformation matrices (which can be time consuming),
       *  Shc::save_transform and Shc::load_transform can be used.
       *
       *  \param[in] transform The transformation matrix
       *  \param[in] matrix_type By specifing the sparsity relations of the matrix, internal computations can be sped up
       *  \param[in] tolerance If a sparse matrix is added, all entries whose absolute values are smaller than \p tolerance can be set to zero
       *  \param[in] l_max The maximal number of bands for which the transformation is computed
       *  \return If successfull, the index of the added transformation is returned. This index is required by Shc::transform(Shpm& shpm, int index).
       *          Otherwise, -1 is returned.
       */
      int add_transform(MatrixReal& transform, e_matrix_type matrix_type = DENSE, float tolerance = 0, int l_max = IGNORE);
      /**
       *  \brief Returns the number of transformations loaded
       *
       *  \return The transformation counter (transformations successfully added via Shc::add_transform and Shc::load_transform)
       */
      int get_transform_size();
      /**
       *  \brief Clears all transformations, the transformation counter is set to zero.
       */
      void clear_transform();
      
      
      
      /**
       *  \brief Creates a weighting function for hemispherial images
       *
       *  Creates a hemispherical weighting function for the \ref compass.
       *  The weighting function only requires the first two bands and is therefore comparably performant.
       *
       *  \return The resulting instance of shc::Shpm
       */
      Shpm create_weighting_hemi();
      /**
       *  \brief Conversion: shpm > shpm
       *
       *  This function can be used to change the \ref hemisphericalcontinuations and maximal number of bands used for an instance of shc::Shpm (Fourier coefficient vector).
       *
       *  \param[in] shpm The original Fourier coefficient vector
       *  \param[in] contin Hemispherical continuation of the result
       *  \param[in] l_max Maximal number of bands of the result
       *  \return The resulting instance of shc::Shpm
       */
      Shpm shpm2shpm(Shpm& shpm, e_contin contin = FULL, int l_max = IGNORE);
      /**
       *  \brief Conversion: surf > shpm
       *
       *  Computes the Fourier coefficient vector of the given surface.
       *
       *  \param[in] surf The surface (panoramic image) for which the Fourier transform is computed
       *  \param[in] contin The used \ref hemisphericalcontinuations 
       *  \param[in] l_max Maximal number of bands used
       *  \return The resulting instance of shc::Shpm
       */
      Shpm surf2shpm(Surf& surf, e_contin contin = FULL, int l_max = IGNORE);
      /**
       *  \brief Conversion: coef > shpm
       *
       *  The class shc::Shpm stores, additionally to the Fourier coefficient vector, additional information internally required by the libShc.
       *  These information are used to efficiently implement various functions, especially if \ref hemisphericalcontinuations are involved.
       *  This function can be used to create an instance of shc::Shpm directly from a Fourier coefficient vector.
       *
       *  \param[in] coef The original Fourier coefficient vector
       *  \param[in] contin The used \ref hemisphericalcontinuations 
       *  \param[in] l_max Maximal number of bands used
       *  \return The resulting instance of shc::Shpm
       */
      Shpm coef2shpm(Coef& coef, e_contin contin = FULL, int l_max = IGNORE);
      /**
       *  \brief Conversion: shpm > coef
       *
       *  Extracts the Fourier coefficient vector from the instance of shc::Shpm.
       *
       *  \param[in] shpm The original instance of shc::Shpm
       *  \return The extracted Fourier coefficient vector
       */
      Coef shpm2coef(Shpm& shpm);
      /**
       *  \brief Conversion: shpm > surf
       *
       *  Computes the inverse Fourier transform
       *
       *  \param[in] shpm The instance of shc::Shpm
       *  \return The inverse Fourier transformed as a surface
       */
      Surf shpm2surf(Shpm& shpm);
      /** 
       *  \brief Conversion: ocamcalib > matrix
       *
       *  In the libShc, it is assumed that panoramic images are used.
       *  However, camera input (e.g. from a fisheye lens) is not a panoramic image.
       *  This function uses the \ref ocamcalib to unwrap camera images directly.
       *  Note that this functionality requires an initialization via Shc::init_ocamcalib.
       *
       *  \param[in] matrix The wrapped camera image
       *  \return The unwrapped panoramic image
       */
      MatrixReal ocamcalib2matrix(MatrixReal& matrix);
      /** \overload
       *  
       *  For higher performance, this function can be used to directly operate on a buffer.
       *       
       *  \param[in] data Pointer to the data buffer in which the wrapped camera image is located
       */
      MatrixReal ocamcalib2matrix(unsigned char* data);
      /** \overload
       *  
       *  For higher performance, this function can be used to directly operate on a buffer.
       *       
       *  \param[in] data Pointer to the data buffer in which the wrapped camera image is located
       *  \param[out] data_out Pointer to the data buffer in which the unwrapped camera image is stored
       */
      void ocamcalib2matrix(unsigned char* data, unsigned char* data_out);
      /** 
       *  \brief Conversion: matrix > surf
       *
       *  Unwraps a camera image using the \ref ocamcalib and creates a surface
       *
       *  \param[in] matrix The wrapped camera image
       *  \return The resulting shc::Surface
       */
      Surf ocamcalib2surf(MatrixReal& matrix);
      /** \overload
       *  
       *  For higher performance, this function can be used to directly operate on a buffer.
       *       
       *  \param[in] data Pointer to the data buffer in which the wrapped camera image is located
       */
      Surf ocamcalib2surf(unsigned char* data);
      /** 
       *  \brief Conversion: ocamcalib > shpm
       *
       *  Unwraps a camera image using the \ref ocamcalib, and performs afterwards the Fourier transform.
       *
       *  \param[in] matrix The wrapped camera image
       *  \param[in] contin \ref hemisphericalcontinuations used for the Fourier transform
       *  \param[in] l_max Maximal number of bands used for the Fourier transform
       *  \return The resulting shc::Shpm
       */
      Shpm ocamcalib2shpm(MatrixReal& matrix, e_contin contin = FULL, int l_max = IGNORE);
      /** \overload
       *  
       *  For higher performance, this function can be used to directly operate on a buffer.
       *
       *  \param[in] data Pointer to a buffer holding the panoramic image information
       *  \param[in] contin \ref hemisphericalcontinuations used for the Fourier transform
       *  \param[in] l_max Maximal number of bands used for the Fourier transform
       */
      Shpm ocamcalib2shpm(unsigned char* data, e_contin contin = FULL, int l_max = IGNORE);
          
      /**
       *  \brief Creates a feature vector containing information about the Fourier transformed function
       *
       *  This functioncan be used to obtain information as for example the amplitude spectrum of a Fourier transformed function.
       *  Note that this representation is sparse, to obtain a specific entry use the corresponding index functions, i.e.
       *  for the bispectrum use Shc::get_feature_index_bs.
       *  Note that for the computation of the bispectrum Clebsch-Gordan matrices are used, see \ref products and Shc::init_bands.
       *  
       *  \param[in] shpm The instance of shc::Shpm for which information is created
       *  \param[in] feature The type of information computed
       *  \return The vector containing the information
       */
      VectorReal get_feature(Shpm& shpm, e_feature feature);
      /**
       *  \brief Compares two panoramic images regarding some measure (e.g. amplitude spectrum)
       *
       *  This function compares two panoramic images by computing for each an information vector (see Shc::get_feature) and computing their difference.
       *  Note that this function is influenced by the settings of Shc::set_feature_norm and Shc::set_feature_normalize
       *  Note that for the computation of the bispectrum Clebsch-Gordan matrices are used, see \ref products and Shc::init_bands.
       *
       *  \param[in] shpm1 The first instance of shc::Shpm
       *  \param[in] shpm2 The second instance of shc::Shpm
       *  \param[in] feature The type of information used for comparison
       *  \return The difference value
       */
      float get_feature_difference(Shpm& shpm1, Shpm& shpm2, e_feature feature);
      /** \overload
       *
       *  This function compares two equally sized sets of shc::Shpm instances at once.
       *
       *  \param[in] shpm1 The first set of shc::Shpm instances
       *  \param[in] shpm2 The second set of shc::Shpm instances
       *  \param[in] feature The type of information used for comparison
       *  \return A vector containing the difference values for each pair
       */
      VectorReal get_feature_difference(VecShpm& shpm1, Shpm& shpm2, e_feature feature);
      /** \brief Computes the index of a specific entry in the information vector
       *
       *  \param[in] shpm The instance of shc::Shpm for which the information vector was computed
       *  \param[in] l Desired band index
       *  \param[in] m Desired band index
       *  \return The index at which the specific entry can be found
       */ 
      int get_feature_index_sh(Shpm& shpm, int l, int m);
      /** \brief Computes the index of a specific entry in the information vector
       *
       *  \param[in] shpm The instance of shc::Shpm for which the information vector was computed
       *  \param[in] l Desired band index
       *  \return The index at which the specific entry can be found
       */ 
      int get_feature_index_as(Shpm& shpm, int l);
      /** \brief Computes the index of a specific entry in the information vector
       *
       *  \param[in] shpm The instance of shc::Shpm for which the information vector was computed
       *  \param[in] l1 Desired band index
       *  \param[in] l2 Desired band index
       *  \param[in] i Desired band index
       *  \return The index at which the specific entry can be found
       */ 
      int get_feature_index_bs(Shpm& shpm, int l1, int l2, int i);
      /** \brief Computes the index of a specific entry in the information vector
       *
       *  \param[in] shpm The instance of shc::Shpm for which the information vector was computed
       *  \param[in] l1 Desired band index
       *  \param[in] l2 Desired band index
       *  \param[in] i Desired band index
       *  \return The index at which the specific entry can be found
       */ 
      int get_feature_index_bs_dense(Shpm& shpm, int l1, int l2, int i);
      
      
      /**
       *  \brief Computes the rotational offset between two panoramic images (visual 3D compass)
       *
       *  Applies the \ref compass approach to determine the rotational offset between two panoramic images directly in the basis of spherical harmonics.
       *  Note that this function is especially influenced by the following aspects:
       *  - \ref rotations
       *  - \ref hemisphericalcontinuations
       *  - \ref tangentdistance
       *
       *  \param[in] shpm1 The <i>current view</i>, which is actively rotated to match the snapshot \p shpm2
       *  \param[in] shpm2 The <i>snapshot</i> with which the current view is compared
       *  \return The rotation which needs to be applied to the current view such that is rotationally aligned with the snapshot
       */
      Xyz compass(Shpm& shpm1, Shpm& shpm2);
      /** \overload
       *
       *  Rotationally aligns two (equally sized) sets of panoramic images.
       *  For each pair (same indices in the vectors), the errors are added up.
       *
       *  \param[in] v_shpm1 A set of <i>current views</i>, which is actively rotated to match the set of snapshots \p v_shpm2
       *  \param[in] v_shpm2 A set of <i>snapshots</i> with which the set of current views is compared
       *  \return The rotation which needs to be applied to the current views such that is rotationally aligned with the snapshots
       */
      Xyz compass(VecShpm& v_shpm1, VecShpm& v_shpm2);
      /** \overload
       *
       *  By using weighting functions (see \ref compass), areas of the panoramic images can be masked out.
       *  The mask should be normalized to the interval [0,1]. 
       *  
       *  \param[in] shpm1 The <i>current view</i>, which is actively rotated to match the snapshot \p shpm2
       *  \param[in] shpm2 The <i>snapshot</i> with which the current view is compared
       *  \param[in] mask1 The mask of the current view
       *  \param[in] mask2 The mask of the snapshot
       *  \return The rotation which needs to be applied to the current view such that is rotationally aligned with the snapshot
       */
      Xyz compass(Shpm& shpm1, Shpm& shpm2, Shpm& mask1, Shpm& mask2);
      /** \overload
       *
       *  Rotationally aligns two sets of panoramic images using masks, the sets of panoramic images and masks are required to have the same size.
       *  For each pair (same indices in the vectors), the errors are added up.
       *  
       *  \param[in] v_shpm1 The set of <i>current views</i>, which is actively rotated to match the set of snapshots \p v_shpm2
       *  \param[in] v_shpm2 The set of <i>snapshots</i> with which the set of current views is compared
       *  \param[in] v_mask1 The masks for each current view
       *  \param[in] v_mask2 The masks for each snapshot
       *  \return The rotation which needs to be applied to the current view such that is rotationally aligned with the snapshot
       */
      Xyz compass(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2);
      /**
       *  /brief Validates the quality of the visual 3D compass estimate
       *  
       *  If the orientations of the current view and snapshot are known, the \ref compass estimate can be validated.
       *
       *  \param[in] rot1 The orientation of the current view
       *  \param[in] rot2 The orientation of the snapshot
       *  \param[in] xyz The compass estimate
       *  \return The rotational offset (smallest angle) between the compass estimate and the real rotational offset
       */
      float compass_evaluate(Xyz rot1, Xyz rot2, Xyz xyz);
      
      /**
       *  \brief Computes a linear combination
       *
       *  Let \f$ v,w \f$ be the Fourier coefficient vectors of \p shpm1 and \p shpm2,
       *  then this function computes the linear combination \f$ u = av + bw \f$.
       *
       * \param[in] a First linear coefficient
       * \param[in] shpm1 First Fourier coefficient vector
       * \param[in] b Second linear coefficient
       * \param[in] shpm2 Second Fourier coefficient vector
       * \return The resulting Fourier coefficient vector
       */
      Shpm linear_combination(float a, Shpm& shpm1, float b, Shpm& shpm2);
      
      /**
       *  \brief Multiplies two panoramic images in the basis of spherical harmonics
       *  
       *  This function computes \ref products between two panoramic images directly in the basis of spherical harmonics.
       *
       *  \param[in] shpm1 The first instance of Shc::Shpm
       *  \param[in] shpm2 The second instance of Shc::Shpm
       *  \return The point-wise product
       */
      Shpm product(Shpm& shpm1, Shpm& shpm2);
      /** \overload
       * 
       *  Computes the poin-wise product between two sets of panoramic images
       *  Note that for the computation of point-wise products, Clebsch-Gordan matrices are used, see Shc::init_bands
       *
       *  \param[in] v_shpm1 The first set of instances of Shc::Shpm
       *  \param[in] v_shpm2 The second set of instances of Shc::Shpm
       *  \return The point-wise products between each pair
       */
      VecShpm product(VecShpm& v_shpm1, VecShpm& v_shpm2);
      /**
       *  \brief Rotates a panoramic image in the basis of spherical harmonics
       *
       *  This function uses the rotations initialized via Shc::init_rotations.
       *  If the initiallized rotations are not sufficient to successfully apply the rotation, a warning is shown.
       *
       *  \param[in] shpm The instance of shc::Shpm to which the rotation is applied
       *  \param[in] xyz The rotation which is applied
       *  \return Returns the rotated panoramic image in the basis of spherical harmonics
       */
      Shpm rotate(Shpm& shpm, Xyz xyz);
      /** \overload
       * 
       *  Computes the rotation for a set of panoramic images
       *
       *  \param[in] v_shpm The set of instances of Shc::Shpm
       *  \param[in] xyz The rotation applied to each panoramic image
       *  \return The rotated panoramic images in the basis of spherical harmonics
       */
      VecShpm rotate(VecShpm& v_shpm, Xyz xyz);
      /**
       *   \brief Translates (Warps) a panoramic image in the basis of spherical harmonics
       *
       *  This function uses the translatons and rotations initialized via Shc::init_rotations and Shc::init_translations.
       *  
       *  See \ref rotations and \ref translations for more details.
       *
       *  \param[in] shpm The instance of shc::Shpm to which the translation is applied
       *  \param[in] translation Target location
       *  \return Returns the translated panoramic image in the basis of spherical harmonics
       */
      Shpm warp(Shpm& shpm, Coor3d translation);
      /** \overload
       * 
       *  Computes the translation for a set of panoramic images
       *
       *  \param[in] v_shpm The set of instances of Shc::Shpm
       *  \param[in] translation Target location
       *  \return The translated panoramic images in the basis of spherical harmonics
       */
      VecShpm warp(VecShpm& v_shpm, Coor3d translation);
      /** 
       *  \brief Applies an arbitrary transformation
       *   
       *  Arbitrary transformations can be created and used to operate on an instance of shc::Shpm.
       *  Each transformation has to be added to the current instance of shc::Shc via Shc::add_transform beforehand.
       *  For more information, see \ref transforms.
       *
       *  \param[in] shpm The current instance of Shc::Shpm
       *  \param[in] index The index of the transformation applied
       *  \return The resulting shc::Shpm
       */
      Shpm transform(Shpm& shpm, int index);
      /** \overload
       *
       *  \param[in] shpm A set of Shc::Shpm instances
       *  \param[in] index The index of the transformation applied
       *  \return The resulting instances of shc::Shpm
       */
      VecShpm transform(VecShpm& shpm, int index);

      /** 
       *  \brief Conversion (OCamCalib): camera coordinate > world coordinate
       *
       *  To use this function, an appropriate calibration file for the currently used camera has to be loaded via Shc::init_ocamcalib.
       *
       *  \param[in] coor2 Coordinate of pixel in the camera image
       *  \result Coordinate of pixel in world coordinates (viewing direction)
       */
      Coor3d ocamcalib_cam2world(Coor2d coor2);
      /** 
       *  \brief Conversion (OCamCalib): world coordinate > camera coordinate
       *
       *  To use this function, an appropriate calibration file for the currently used camera has to be loaded via Shc::init_ocamcalib.
       *
       *  \param[in] coor3 Coordinate of pixel in world coordinates (viewing direction)
       *  \result Coordinate of pixel in the camera image
       */
      Coor2d ocamcalib_world2cam(Coor3d coor3);
      
      /**
       *  \brief Saves all currently loaded transformations to a file
       *
       *  Saves all transformations successfully added via Shc::add_transform and Shc::load_transform to the specified file.
       *
       *  \param[in] file Filename
       *  \return Returns if the transformation matrices could successfully be saved to the specified file
       */
      bool save_transform(std::string file);
      /**
       *  \brief Loads transformation matrices from a file
       *
       *  Loads transformations from the specified file and adds them to the current instance of shc::Shc.
       *
       *  \param[in] file Filename
       *  \return Returns if the transformation matrices could successfully be loaded from the specified file
       */ 
      bool load_transform(std::string file);
      bool save_tangent_distance(std::string file);
      bool load_tangent_distance(std::string file);

      /** 
       *  \brief Returns the number of bands used by this instance
       *
       *  Returns the number of bands set via Shc::init_bands.
       *
       *  \param[out] n_bands The number of bands
       *  \param[out] n_bands_CG The number of bands for Clebsch-Gordan matrices
       */
      void get_bands(int& n_bands, int& n_bands_CG);
      
      
  };

}



#endif











