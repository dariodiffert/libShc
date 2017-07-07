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

void Shc::filter_bilinear_indices(float& y, float& x, float h, float w, int& y0, int& x0, int& y1, int& x1) {
  
  x0 = floor(x); x1 = x0+1;
  y0 = floor(y); y1 = y0+1;
    
  if (x0 < 0) {
    x0 = 0;
  } else if (x0 > w-1) {
    x0 = w-1;
  }
  if (x1 < 0) {
    x1 = 0;
  } else if (x1 > w-1) {
    x1 = w-1;
  }
  if (y0 < 0) {
    y0 = 0;
  } else if (y0 > h-1) {
    y0 = h-1;
  }
  if (y1 < 0) {
    y1 = 0;
  } else if (y1 > h-1) {
    y1 = h-1;
  }
  
  x = x - x0;
  y = y - y0;
  
}
void Shc::filter_bilinear_getValue(MatrixReal& matrix, int y0, int x0, int y1, int x1, float& f00, float& f10, float& f01, float& f11) {
  
  f00 = matrix(y0,x0);
  f10 = matrix(y1,x0);
  f01 = matrix(y0,x1);
  f11 = matrix(y1,x1);
  
}
void Shc::filter_bilinear_getValue(unsigned char* data, int width, int y0, int x0, int y1, int x1, float& f00, float& f10, float& f01, float& f11) {
    
  f00 = (float) data[y0*width+x0];
  f10 = (float) data[y1*width+x0];
  f01 = (float) data[y0*width+x1];
  f11 = (float) data[y1*width+x1];
  
}
float Shc::filter_bilinear_calc(float f00, float f10, float f01, float f11, float y, float x) {
    
  float result = f00*(1-y)*(1-x) + f10*y*(1-x)+f01*(1-y)*x+f11*y*x;
  
  return result;
  
}
float Shc::filter_bilinear_mat(MatrixReal& matrix, int h, int w, float y, float x) {
  
  int x0, x1, y0, y1;
  float f00, f10, f01, f11;
  
  filter_bilinear_indices(y, x, h, w, y0, x0, y1, x1);
  filter_bilinear_getValue(matrix, y0, x0, y1, x1, f00, f10, f01, f11);
  float result = filter_bilinear_calc(f00, f10, f01, f11, y, x);
  
  return result;
  
}
float Shc::filter_bilinear_raw(unsigned char* data, int h, int w, float y, float x) {
  
  int x0, x1, y0, y1;
  float f00, f10, f01, f11;
  
  filter_bilinear_indices(y, x, h, w, y0, x0, y1, x1);
  filter_bilinear_getValue(data, w, y0, x0, y1, x1, f00, f10, f01, f11);
  float result = filter_bilinear_calc(f00, f10, f01, f11, y, x)  / 255.0f;
  
  return result;
  
}

void Shc::filter_nearest_indices(float y, float x, float h, float w, int& y0, int& x0) {
  
  x0 = x + 0.5;
  y0 = y + 0.5;
  
  if (x0 >= w-1) {
    x0 = w-1;
  }
  if (y0 >= h-1) {
    y0 = h-1;
  }
  if (x0 < 0) {
    x0 = 0;
  }
  if (y0 < 0) {
    y0 = 0;
  }
  
}
float Shc::filter_nearest_mat(MatrixReal& matrix, int h, int w, float y, float x) {
  
  int x0, y0;
  filter_nearest_indices(y, x, h, w, y0, x0);
  float result = matrix(y0,x0);
  
  return result;
  
}
float Shc::filter_nearest_raw(unsigned char* data, int h, int w, float y, float x) {
  
  int x0, y0;
  filter_nearest_indices(y, x, h, w, y0, x0);
  float result = data[y0*w+x0] / 255.0f;
  
  return result;
  
}

