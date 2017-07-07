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

VectorReal Shc::stepspace(float resolution, float max_val, bool negative) {
  
  vector<float> v;
  
  v.push_back(0);
  
  float x = resolution;
  while (true) {
    v.push_back(x);
    if (negative) {
      v.push_back(-x);
    }
    if (x >= max_val-1e-5) {
      break;
    } else {
      x+=resolution;
    }
  }
  
  int n = v.size();
  VectorReal result(n);
  for (int i=0; i<n; i++) {
    result(i) = v[i];
  }
  
  return result;
  
}
void Shc::equalize_vector_bigger(Shpm& shpm1, Shpm shpm2, Coef& v1, Coef& v2) {
  
  // !!! remove later ???
  if (shpm1.coef.size() != l2sh(shpm1.l_max, shpm1.contin)) {
    print_warning("equalize_size", "Input vector has invalid size! BUG");
  }
  
  // !!! remove later ???
  if (shpm2.coef.size() != l2sh(shpm2.l_max, shpm2.contin)) {
    print_warning("equalize_size", "Input vector has invalid size! BUG");
  }
  
  v1 = shpm1.coef;
  v2 = shpm2.coef;
  
  int n1 = v1.size();
  int n2 = v2.size();
  
  if (n1 < n2) {
    v1.conservativeResizeLike(VectorReal::Zero(n2));
  } else if (n2 < n1) {
    v2.conservativeResizeLike(VectorReal::Zero(n1));
  }
  
}
float Shc::min3(float i1, float i2, float i3) {
  return min(min(i1,i2),i3);
}
float Shc::max3(float i1, float i2, float i3) {
  return max(max(i1,i2),i3);
}
Shpm Shc::linear_combination(float a, Shpm& shpm1, float b, Shpm& shpm2) {
  
  Shpm result;
  Shpm* ref1;
  Shpm* ref2;
  
  Shpm shpm1_full, shpm2_full;
  
  if (trifold(shpm1.contin) == trifold(shpm2.contin)) {
    ref1 = &shpm1;
    ref2 = &shpm2;
    result.contin = shpm1.contin;
  } else {
    shpm1_full = contin_convert(shpm1, FULL);
    shpm2_full = contin_convert(shpm2, FULL);
    ref1 = &shpm1_full;
    ref2 = &shpm2_full;
    result.contin = FULL; 
  }
  result.l_max = max(shpm1.l_max, shpm2.l_max);
  
  int n1 = ref1->coef.size();
  int n2 = ref2->coef.size();
  
  if (n1 < n2) {
    Coef c1(n2);
    c1.segment(0, n1) = ref1->coef;
    c1.segment(n1, n2-n1) = VectorReal::Zero(n2-n1);
    result.coef = a*c1+b*ref2->coef;
  } else if (n2 < n1) {
    Coef c2(n1);
    c2.segment(0, n2) = ref2->coef;
    c2.segment(n2, n1-n2) = VectorReal::Zero(n1-n2);
    result.coef = a*ref1->coef+b*c2;
  } else {
    result.coef = a*ref1->coef+b*ref2->coef;
  }

  return result;
  
}

VectorReal Shc::intern_calcHistogram(MatrixReal& matrix, float theta_min, float theta_max, bool cumulative, bool panoramic) {
  
  int m = matrix.rows();
  int n = matrix.cols();
  VectorReal hist(256);
  hist.setZero();
  
  for (int y=0; y<m; y++) {
    for (int x=0; x<n; x++) {
      int h = round(matrix(y,x) * 255.0);
      float t = (float)(y)/(float)(m) * (theta_max - theta_min) + theta_min;
      float w;
      if (panoramic == true) {
        w = sin(t);
      } else {
        w = 1;
      }
      if (h<0) {
        h=0;
      }
      if (h>255) {
        h=255;
      }
      hist(h) += w;
    }
  }

  if (cumulative == true) {
    for (int i=1; i<256; i++) {
      hist(i) += hist(i-1);
    }
  }
  
  return hist;
  
}
float Shc::intern_histogramEqualization_calcEntry(VectorReal& hist, float val) {
  
  int c = round(val * 255.0);
  
  if (c<0) {
    c = 0;
  }
  if (c>255) {
    c = 255;
  }

  float result = (hist(c) - hist(0)) / (hist(255) - hist(0));
  return result;
  
}
float Shc::intern_otsu(VectorReal& hist) {
  
  float total = 0;
  for (int i=0; i<256; i++) {
    total += hist[i];
  }
  
  float t       = 0;
  float sumB    = 0;
  float wB      = 0;
  float maximum = 0;
  float sum1    = 0;
  for (int i=0; i<256; i++) {
    sum1 += i*hist[i];
  }
  
  for (int i=0; i<256; i++) {
    wB = wB + hist[i];
    if (wB == 0) {continue;}
    float wF = total - wB;
    if (wF == 0) {break;}
    sumB = sumB +  i * hist[i];
    float mB = sumB / wB;
    float mF = (sum1 - sumB) / wF;
    float between = wB * wF * (mB - mF) * (mB - mF);
    if (between >= maximum) {
      t = i;
      maximum = between;
    }
  }
  
  return t / 255.0f;
  
}



