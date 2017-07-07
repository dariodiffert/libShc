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

// Image class wrapper functions. In this case the image container is a "cv::Mat" from the openCV library.
// Adapt the four functions below and the typedef in Shx.h (Image) to use a custom image library.

void Shx::show(Image& image, string title, bool waitKey) {
  // show image in a window with given title.
  // Optionally wait for a key press to continue.
  
  imshow(title, image);
  if (waitKey) {
    print("show: press any key to continue...\n");
    cv::waitKey(0);
  } else {
     cv::waitKey(1);
  }
  
}
bool Shx::read_image(Image& image, string file) {
  // Read an image from a file.

  image = cv::imread(file, CV_LOAD_IMAGE_GRAYSCALE);
  if(!image.data ) {
    return false;
  } else {
    return true;
  }
  
}
bool Shx::write_image(Image& image, string file) {
  // Write an image to a file.

  return cv::imwrite(file, image);
  
}
MatrixReal Shx::image2matrix(Image &in) {
  // Converts the Image object into a MatrixReal (typedef of MatrixXf) from the Eigen library.
  // This is done (for simplicity) by reading the pixels one-by-one.
  // Note that the input gray values are divided by 255.0 to ensure that internally all values range between [0,1].
  // For time critical applications access surfaces directly (see manual).
  
  int h = in.rows;
  int w = in.cols;
  
  MatrixReal out(h, w);
  
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      // it is advantageous to normalize the color values to the interval [0,1] to make multiplications (e.g. with binary masks) more intuitive.
      out(y,x) = in.at<uchar>(y,x) / 255.0;
    }
  }
  
  return out;
  
}
Image Shx::matrix2image(MatrixReal &in) {
  // the internal representation via an Eigen matrix is converted back to an openCV image (compare image2matrix).
  // Note the boundary checks and normalization to ensure that all values are in [0,255].
  
  int h = in.rows();
  int w = in.cols();
  
  Image out(h, w, CV_8UC1);
  
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      
      // undo the normalization from image2matrix()
      float val = in(y,x);

      uchar ret;
      if (val <= 0) {
        ret = 0;
      } else if (val >= 1) {
        ret = 255;
      } else {
        ret = val * 255.0;
      }
      
      out.at<uchar>(y,x) = ret;
        
    }
  }
  
  return out;
  
}

 
