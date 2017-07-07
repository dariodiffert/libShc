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

#ifndef AGTI_DD_SHC_SHX
#define AGTI_DD_SHC_SHX

#include "Shc.h"

// --- EDIT HERE --- custom image library (Shc itself does not support any image handling)
#include <opencv2/opencv.hpp>
// --- EDIT HERE ---

namespace shc {

  // --- EDIT HERE --- Adapt to use custom image container
  /**
   *  The image container used to load/save/display images.
   *  By default, we use the OpenCV (\ref installation) for these tasks,
   *  however by adjusting Shx_wrapper.cpp, arbitrary image libraries can be used instead.
   */
  typedef cv::Mat Image; 
  // --- EDIT HERE ---
   
   /**
   *   \brief An extension of shc::Shc for loading/writing images from the disk or showing images.
   * 
   * This class provides several functions to load and write images to the disk as well as images to show images on the screen.
   * To implement this functionality, the OpenCV (\ref installation) is used.
   * For an overview of all classes and their relations, see \ref classlayout.
   */ 
  class Shx : public Shc {
    
    private:
      
      // wrapper functions
      void show(Image& image, std::string title = "", bool waitKey = true);
      bool read_image(Image& image, std::string file);
      bool write_image(Image& image, std::string file);
      Image matrix2image(MatrixReal &in);
      MatrixReal image2matrix(Image &in);  
      Image load_image(std::string file);
      bool save_image(Image& image, std::string file);
      
    public:
      
      /** 
       *  \brief Loads a matrix from an panoramic image file
       *  \param[in] file Path to file on disk
       *  \return Returns the loaded matrix
      */  
      MatrixReal load_matrix(std::string file);
      
      /** 
       *  \brief Saves a matrix as a panoramic image file
       *  \param[in] matrix Matrix to save
       *  \param[in] file Path to file on disk
       *  \return Returns if the matrix could successfully be written to the file
      */  
      bool save_matrix(MatrixReal& matrix, std::string file);
      
      /** 
       *  \brief Loads a shc::Shpm instance from a panoramic image file
       *  \param[in] file Path to file on disk
       *  \param[in] contin Hemispherical continuation used
       *  \param[in] l_max Maximal number of bands used for Fourier transform
       *  \return Returns the loaded instance of shc::Shpm.
      */  
      Shpm load_shpm(std::string file, e_contin contin = FULL, int l_max = IGNORE);
      
      /** 
       *  \brief Saves a shc::Shpm instance as a panoramic image file 
       *  \param[in] shpm Instance to save
       *  \param[in] file Path to file on disk
       *  \return Returns if the shc::Shpm instance could successfully be written to the file
      */  
      bool save_shpm(Shpm& shpm, std::string file);

      /** 
       *  \brief Loads a shc::Shpm instance from a camera image file. The camera image is unwrapped using the OCamCalib
       *  \param[in] file Path to file on disk
       *  \param[in] contin Hemispherical continuation used
       *  \param[in] l_max Maximal number of bands used for Fourier transform
       *  \return Returns the loaded instance of shc::Shpm.
       */  
      Shpm load_shpm_ocamcalib(std::string file, e_contin contin = FULL, int l_max = IGNORE);
      
      /** 
       *  \brief Plots the matrix as a panoramic image in a window using OpenCV.
       *  \param[in] matrix The matrix which should be plotted
       *  \param[in] title An optional parameter to set the window title. By using the same title again, the window is redrawn with the new matrix.
       *  \param[in] waitKey If set, the execution is stopped until the user presses a key. The window needs to have focus to register the key press.
      */  
      void show(MatrixReal& matrix, std::string title = "", bool waitKey = true);     
      
      /** \overload 
       *  \brief Plots the inverse Fourier transform of the shc::Shpm instance
      */
      void show(Shpm& shpm, std::string title = "", bool waitKey = true);
      
      /** \overload
       *  \brief Plots a list of matrices
      */
      void show(std::initializer_list<MatrixReal> matrix, std::string title = "", bool waitKey = true);
      
      /** \overload
       *  \brief Plots a list of shc::Shpm instances
      */
      void show(std::initializer_list<Shpm> shpm, std::string title = "", bool waitKey = true);
      
      /** \overload
       *  \brief Plots a vector containing matrices
      */
      void show(std::vector<MatrixReal>& matrix, std::string title = "", bool waitKey = true);
      
      /** \overload
       *  \brief Plots a vector containing shc::Shpm instances
      */
      void show(std::vector<Shpm>& shpm, std::string title = "", bool waitKey = true);
      
  };

}

#endif


