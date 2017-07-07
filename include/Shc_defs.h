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

#ifndef AGTI_DD_SHCDEFS
#define AGTI_DD_SHCDEFS

namespace shc {

  // ----------------------------------- ENUMS & CONSTANTS -----------------------------------
  
  static const int IGNORE = -1;
  
  
  /** 
   *  Determines an axis as used in the \ref coordinatesystem.
   */  
  enum e_axis {
    AXIS_X, ///< The X-axis
    AXIS_Y, ///< The Y-axis
    AXIS_Z  ///< The Z-axis
  };
  /** 
   *  The basis of real spherial harmonics is defined for the complete sphere, therefore only functions which are well-defined on the complete sphere can be Fourier transformed.
   * However, in practical applications (e.g. hemispherical cameras) a function might be defined for a hemisphere only.
   * One possibility is to fill the missing information by random or constant data (see \ref hemisphericalcontinuations).
   * Due to the symmetry of the spherical harmonics this can also be implemented by mirroring the upper hemisphere on the lower.
   * Depending on the mirroring approach this can be implemented in a highly efficient way.
   * Using HEMI_RM and HEMI_RMN reduces the computation times; the computation of translations and the bispectrum speeds up by factor 4,
   * most other calculations (e.g. rotations, amplitude spectrum, etc) speed up by factor 2.
   */ 
  enum e_contin {
    FULL,    ///< The complete sphere is transformed.
    HEMI_M,  ///< The upper hemisphere is mirrored down.
    HEMI_MN, ///< The upper hemisphere is mirrored down and negated.
    HEMI_RM, ///< The upper hemisphere is mirrored down and rotated by \f$ 180^\circ \f$ degree.
    HEMI_RMN ///< The upper hemisphere is mirrored down, rotated by \f$ 180^\circ \f$, and negated.
  };
  /** 
   * The function Shc::get_feature_difference calculates the difference between two vectors \f$ v_1, v_2 \f$. The used norm can be one of the following.
   */
  enum e_feature_norm {
    NORM_1,       ///< Absolute difference \f$ |v_1 - v_2| \f$
    NORM_2,       ///< Euclidean norm \f$ \|v_1 - v_2\| \f$
    NORM_INFINITY ///< Maximum norm \f$ \|v_1 - v_2\|_\infty \f$
  }; 
  /** 
   * A total of three different features can be obtained from an arbitrary Fourier transformed function.
   * These features can be used to get information over a Fourier transformed function \f$ f \f$ (Shc::get_feature) or compare two functions \f$ f,g \f$ (Shc::get_feature_difference).
   */
  enum e_feature {
    ISE,     ///< The <i>Integral Squared Error</i> between two functions is defined as \f$ \int_{S^2} (f(s)-g(s))^2 ds = \|a_f-a_g\|^2 \f$.
    AS,      ///< The <i>Amplitude Spectrum</i> is calculated.
    BS,      ///< The <i>Bispectrum</i> is calculated.
    BS_DENSE ///< The <i>Bispectrum</i> is calculated, however only for all pairs of bands with \f$ |l_1 - l_2| \le 1 \f$.
  };
  /** 
   * \ref translations can be simulated directly in the basis of real spherical harmonics. Depending on the approach, a translation of a function \f$ f \f$ can be interepreted in different ways. The following two are currently implemented:
   */
  enum e_trans_type {
    VISUAL, ///< Let the function \f$ f \f$ represent a visual scene projected on a sphere. Then a translation of the center is equivalent to a translation of the viewpoint. Using the <b>equal distance assumption</b> the change of the visual scene can be predicted.
    DENSITY ///< For a function \f$ f \f$  representing a density (e.g. a \ref pointclouds) a translation does not only change the center, but also the appearance of the density itself (assuming that the density is evenly distributed): If the center moves towards the sphere, the observed density in the direction of motion will decrease, while it increases in the opposite direction. The increase and decrease depends on both the translation distance and the observed angle on the sphere.
  };
  /** 
   * The shc::Shc and shc::Shx classes produce automaticly output and print it to the console. This output can be limited if necessary:
   */
  enum e_print_level{
    NONE,     ///< All information and warnings are printed.
    WARNINGS, ///< Only warnings are printed.
    ALL       ///< Nothing is printed, note that ignoring warnings may lead to undefined behaviour! User calls to Shc::print are always shown.
  };
  /**  
   * Whenever image input is used (e.g. Shx::load_shpm) it can be chosen between the following two filters to interpolate floating point image coordinates:
   */
  enum e_filter {
    NEAREST, ///< Takes the nearest/closest pixel.
    BILINEAR ///< Uses bilinear interpolation.
  };
  /** 
   * \ref transforms are expressed by arbitrary transformation matrices \f$ \mathbf{M} \f$.
   * Internal calculations can be sped up if the form of the matrix is known:
   */ 
  enum e_matrix_type {
    DENSE,   ///< All/most entries of \f$ \mathbf{M} \f$ are non-zero.
    SPARSE,  ///< Most entries of \f$ \mathbf{M} \f$ are zero.
    BANDWISE ///< The matrix \f$ \mathbf{M} = \bigoplus_i \mathbf{M}_i \f$ consists of band-wise block matrices, where the i-th block has dimension \f$ 2i+1 \times 2i+1 \f$ (\ref rotations).
  };
  /** 
   * Determines the modus of the \ref tangentdistance.
   */
  enum e_tangent_distance_mode {
    OFF,      ///< Deactivates the tangent distance.
    ONESIDED, ///< The one-sided tangent distance is calculated.
    TWOSIDED  ///< The two-sided tangent distance is calculated.
  };
  /** 
   * Used to apply a predefined convolution via Shc::convolution.
   * The implemented convolution types are:
   * <table>
   * <tr><th> e_convolution_type <th> kernel <th> type
   * <tr><td> SOBEL         <td> \f$ \begin{pmatrix} 1 & 0 & -1 \\ 2 & 0 & -2 \\ 1 & 0 & -1 \end{pmatrix} \f$  & \f$ \begin{pmatrix} 1 & 2 & 1 \\ 0 & 0 & 0 \\ -1 & -2 & -1 \end{pmatrix} \f$ <td> Edge detection
   * <tr><td> PREWITT       <td> \f$ \begin{pmatrix} 1 & 0 & -1 \\ 1 & 0 & -1 \\ 1 & 0 & -1 \end{pmatrix} \f$ & \f$ \begin{pmatrix} 1 & 1 & 1 \\ 0 & 0 & 0 \\ -1 & -1 & -1 \end{pmatrix} \f$  <td> Edge detection
   * <tr><td> EDGE          <td> \f$ \begin{pmatrix} 0 & 1 & 0 \\ 1 & -4 & 1 \\ 0 & 1 & 0 \end{pmatrix} \f$            <td> Edge detection
   * <tr><td> LAPLACE       <td> \f$ \begin{pmatrix} 1 & 2 & 1 \\ 2 & -12 & 2 \\ 1 & 2 & 1 \end{pmatrix} \f$           <td> Edge detection
   * <tr><td> BLUR_BOX      <td> \f$ \frac{1}{9}\begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \end{pmatrix} \f$  <td> Blur
   * <tr><td> BLUR_GAUSSIAN <td> \f$ \frac{1}{16}\begin{pmatrix} 1 & 2 & 1 \\ 2 & 4 & 2 \\ 1 & 2 & 1 \end{pmatrix} \f$ <td> Blur
   * </table>
   * The result is then given by
   * \f$ \mathbf{I'} = |\mathbf{I} * \mathbf{K}| \f$ or \f$ \mathbf{I'} = |\mathbf{I} * \mathbf{K_1}| + |\mathbf{I} * \mathbf{K_2}| \f$ 
   * for single and dual kernels, respectively, where the absolut operation is applied elementwise.
   */
  enum e_convolution_type {
    SOBEL,         ///< See table above
    BLUR_BOX,      ///< See table above
    BLUR_GAUSSIAN, ///< See table above
    PREWITT,       ///< See table above
    EDGE,          ///< See table above
    LAPLACE        ///< See table above
  };
  /** 
   * If full-spherical panoramic images are not available, missing parts can be filled-in using noise.
   * The type of noise is precalculated during the initialization phase using Shc::init_noise.
   * There are four different types of noise which can be used, which differ by their amplitude spectrum.
   */
  enum e_noise {
    CONSTANT, ///< The amplitude spectrum is constant.
    NATURAL,  ///< The amplitude spectrum follows a natural distribution (see \ref surface).
    ZERO,     ///< The amplitude spectrum is zero (black).
    CUSTOM    ///< Noise is created which has the amplitude spectrum passed to Shc::init_noise.
  };
  
  // ----------------------------------- CLASSES & STRUCTS -----------------------------------
  
  class Shc;
  
  /**
   *   @brief Stores a Fourier coefficient vector and additional information.
   * 
   * This class is used to store all the Fourier coefficient vector of a function, the maximal number of bands (see \ref fouriertransform), and the used \ref hemisphericalcontinuations.
   * These additional information are needed for efficient computations via the class shc::Shc.
   * Note that the Fourier coefficient vector cannot be directly accessed, since it can be a sparse representation.
   * Instead use Shc::shpm2coef to obtain the Fourier coefficient vector.
   * 
   */ 
  class Shpm {
    friend class Shc;
    public:
      /// @brief Default constructor.
      Shpm() {contin = FULL; l_max = 0;};
      /**
       *  @brief Returns if the instance is initialized, i.e. a Fourier coefficient vector is stored.
       * 
       *  \return If true, the instance is initialized.
       */ 
      bool initialized() {return l_max > 0 ? true : false;};
      /**
       *  @brief Returns the maximal band \f$ L \f$ of the stored Fourier coefficient vector.
       * 
       *  \return Maximal band index \f$ L \f$.
       */ 
      int get_max_band() {return l_max;}
      /**
       *  @brief Returns the \ref hemisphericalcontinuations used.
       * 
       *  \return hemispherical continuation mode used. If no hemispherical continuation is used, it returns <i>FULL</i>.
       */ 
      e_contin get_contin() {return contin;}
    private:
      Eigen::VectorXf coef;
      int l_max;
      e_contin contin;
  };
  
  /**
   *   @brief Angle/axis representing of rotations
   * 
   * This class can be used to represent rotations via its rotation angle and rotation axis.
   * Make sure to use the correct \ref coordinatesystem.
   */ 
  class Axr {
    friend class Shc;
    public:
      
      /**
       *   @brief Rotation angle in radians.
       */ 
      float angle;
      /**
       *  @brief Rotation axis.
       *  The rotation axis is represented by a vector.
       *  The vector length is not considered for conversions via Shc::axr2r or similar, i.e. it is not required to be a unit vector.
       */ 
      Eigen::Vector3f axis;
      
      /// @brief Default constructor.
      Axr() : angle(0), axis(1,0,0) {};
      /**
       *  @brief Constructor used to represent a rotation via its rotation angle and axis.
       * 
       *  \param[in] angle Rotation angle in radians.
       *  \param[in] axis Rotation axis, see \ref coordinatesystem.
       *  \param[in] deg If set, \p angle is interpreted in degrees.
       */ 
      Axr(float angle, Eigen::Vector3f axis, bool deg = false) {
        if (deg) {
          this->angle = angle*M_PI/180.0;
        } else {
          this->angle = angle;
        }
        this->axis = axis;
      }
    
  };
  
  /**
   *   @brief Representation of a rotation via the X"Y'Z Tait-Bryan convention.
   * 
   *  This class is used to represent rotations via the X"Y'Z Tait-Bryan convention,
   *  that is a rotation around the Z-axis, followed by a rotation around the Y-axis, and finally a rotation around the X-axis.
   *  Make sure to use the correct \ref coordinatesystem.
   */ 
  class Xyz {
    friend class Shc;
    public:
      
      /// @brief Rotation angle around the X-axis in radians
      float x;
      /// @brief Rotation angle around the Y-axis in radians
      float y;
      /// @brief Rotation angle around the Z-axis in radians
      float z;
    
      /// @brief Default constructor.
      Xyz() : x(0), y(0), z(0) {};
      /**
       *  @brief Constructor used to represent a rotation via its X"Y'Z Tait-Bryan convention.
       *  
       *  \param[in] x Rotation around the X-axis in radians. The X-axis rotation is applied last.
       *  \param[in] y Rotation around the Y-axis in radians. The Y-axis rotation is applied second.
       *  \param[in] z Rotation around the Z-axis in radians. The Z-axis rotation is applied first.
       *  \param[in] deg If set, \p x, \p y, and \p z are interpreted in degrees.
       */
      Xyz(float x, float y, float z, bool deg = false) {
        if (deg) {
          this->x = x*M_PI/180.0; this->y = y*M_PI/180.0; this->z = z*M_PI/180.0;
        } else {
          this->x = x; this->y = y; this->z = z;
        }
      };
  };
   
  
  /** 
   * A vector of integer values (derived from Eigen, see \ref installation).
   */  
  typedef Eigen::VectorXi                            VectorInt;
  /** 
   * A vector of float values (derived from Eigen, see \ref installation).
   */  
  typedef Eigen::VectorXf                            VectorReal;
  /** 
   * A vector of float values (derived from Eigen, see \ref installation).
   * Used to store Fourier coefficient vectors.
   */  
  typedef Eigen::VectorXf                            Coef;
  /** 
   * A vector of float values (derived from Eigen, see \ref installation).
   * Used to store surface information, see \ref surface.
   */  
  typedef Eigen::VectorXf                            Surf;
  /** 
   * A vector of float values with three entries (derived from Eigen, see \ref installation).
   * Used to store a 3D coordinate.
   */  
  typedef Eigen::Vector3f                            Coor3d;
  /** 
   * A vector of float values with two entries (derived from Eigen, see \ref installation).
   * Used to store a 2D coordinate.
   */  
  typedef Eigen::Vector2f                            Coor2d;
  /** 
   * A matrix of float values (derived from Eigen, see \ref installation).
   * Primarily used to represent and manipulate panoramic images.
   */  
  typedef Eigen::MatrixXf                            MatrixReal;
  /** 
   * A 3x3 matrix of float values (derived from Eigen, see \ref installation).
   * Used to represent rotations, see \ref rotationparameters.
   */  
  typedef Eigen::Matrix3f                            MatrixRotation;
  /** 
   * A vector of shc::Shpm instances.
   * Primarily used to handle multiple slices, see \ref translations.
   */  
  typedef std::vector<Shpm>                          VecShpm;
  /** 
   * A vector of 3D coordinates (see above).
   * Used to represent \ref pointclouds.
   */  
  typedef std::vector<Coor3d>                       Pointcloud;
    
  /// \cond
  typedef std::vector<Eigen::VectorXf> VecCoef;
  typedef std::vector<Eigen::VectorXf> VecSurf;
  typedef std::vector<Eigen::VectorXf> CoefBandwise;
  typedef std::vector<Eigen::MatrixXf> VecMatrixReal;
  typedef std::vector<Eigen::VectorXi> VecVectorInt;
  typedef Eigen::SparseMatrix<float> MatrixRealSparse;
  typedef Eigen::SparseMatrix<std::complex<float> > MatrixComplexSparse;
  /// \endcond
  
  
}
  
#endif











