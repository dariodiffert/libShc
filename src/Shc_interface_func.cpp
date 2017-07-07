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

Pointcloud Shc::load_pointcloud_auto(string file_pointcloud) {
  
  if (!check_init()) {return Pointcloud();}
  
  Pointcloud pointcloud = load_pointcloud(file_pointcloud);
  
  int n = pointcloud.size();
  Coor3d center(0,0,0);
  for (int i=0; i<n; i++) {
    center += pointcloud[i];
  }
  center = center / n;
  
  pointcloud = warp(pointcloud, center);
  
  float curNorm;
  float maxNorm = 0;
  for (int i=0; i<n; i++) {
    curNorm = pointcloud[i].norm();
    if (curNorm > maxNorm) {
      maxNorm = curNorm;
    }
  }
  
  float maxRadius = lut_translate.slices(lut_translate.n_slices-1);
  
  if (maxNorm > 0) {
    pointcloud = warp(pointcloud, Coor3d(0,0,0), maxRadius / maxNorm);
  }
  
  return pointcloud;
  
}
Pointcloud Shc::load_pointcloud(string file_pointcloud) {

  if (!check_init()) {return Pointcloud();}
  
  return read_file_vec3(file_pointcloud);

}
Pointcloud Shc::rotate(Pointcloud& pointcloud, Xyz rotation) {
  
  if (!check_init()) {return Pointcloud();}
  
  int n = pointcloud.size();
  Pointcloud res(n);
  
  MatrixRotation R = xyz2r(rotation);
  for (int i=0; i<n; i++) {
    res[i] = R * pointcloud[i];
  }
  
  return res;
  
}
Pointcloud Shc::warp(Pointcloud& pointcloud, Coor3d translation, float scale) {
  
  if (!check_init()) {return Pointcloud();}
  
  int n = pointcloud.size();
  Pointcloud res(n);
  
  for (int i=0; i<n; i++) {
    res[i] = (pointcloud[i] - translation) * scale;
  }
  
  return res;
  
}

Xyz Shc::compass(Shpm& shpm1, Shpm& shpm2) {
  
  VecShpm s1;
  VecShpm s2;
  VecShpm m1;
  VecShpm m2;
  
  s1.push_back(shpm1);
  s2.push_back(shpm2);
  
  return compass(s1, s2, m1, m2);
  
}
Xyz Shc::compass(Shpm& shpm1, Shpm& shpm2, Shpm& mask1, Shpm& mask2) {
  
  VecShpm s1;
  VecShpm s2;
  VecShpm m1;
  VecShpm m2;
  
  s1.push_back(shpm1);
  s2.push_back(shpm2);
  m1.push_back(mask1);
  m2.push_back(mask2);
  
  return compass(s1, s2, m1, m2);
  
}
float Shc::compass_evaluate(Xyz rot1, Xyz rot2, Xyz xyz) {
  
  float result = angular_difference(product(xyz, rot1), rot2);
  
  return result;
  
}

Shpm Shc::product(Shpm& shpm1, Shpm& shpm2) {

  VecShpm v_shpm1, v_shpm2;
  v_shpm1.push_back(shpm1);
  v_shpm2.push_back(shpm2);
  
  VecShpm v_result = product(v_shpm1, v_shpm2);
  
  return v_result[0];
  
}
Shpm Shc::rotate(Shpm& shpm, Xyz xyz) {
  
  VecShpm v_shpm;
  v_shpm.push_back(shpm);
  
  VecShpm v_result = rotate(v_shpm, xyz);
  
  return v_result[0];
  
}
Shpm Shc::warp(Shpm& shpm, Coor3d translation) {
  
  VecShpm v_shpm;
  v_shpm.push_back(shpm);
  
  VecShpm v_result = warp(v_shpm, translation);
  
  return v_result[0];
  
}
Shpm Shc::transform(Shpm& shpm, int index) {
  
  VecShpm v_shpm;
  v_shpm.push_back(shpm);
  
  VecShpm result = transform(v_shpm, index);

  return result[0];
  
}

VecShpm Shc::product(VecShpm& v_shpm1, VecShpm& v_shpm2) {

  if (!check_init() || !check_init(v_shpm1) || !check_init(v_shpm2)) {return v_shpm1;}
  
  if (equal(v_shpm1.size(), v_shpm2.size()) == false) {
    print_warning("product", "Input vectors have different sizes");
    return v_shpm1;
  }
   
  VecShpm v_result = intern_product(v_shpm1, v_shpm2);
  
  return v_result;
  
}
VecShpm Shc::rotate(VecShpm& v_shpm, Xyz xyz) {
  
  if (!check_init() || !check_init(v_shpm)) {return v_shpm;}
 
  VecShpm v_res = intern_rotate(v_shpm, xyz);
  
  return v_res;
  
}
VecShpm Shc::warp(VecShpm& v_shpm, Coor3d translation) {
  
  if (!check_init() || !check_init(v_shpm)) {return v_shpm;}
  
  int n = v_shpm.size();
  if (n != lut_translate.n_slices) {
    print_warning("warp", "Size of input vectors has to equal the number of slices.");
    return v_shpm;
  }
  
  VecShpm v_result = intern_warp(v_shpm, translation);
  
  return v_result;
  
}
VecShpm Shc::transform(VecShpm& v_shpm, int index) {
  
  if (!check_init()) {return v_shpm;}
  if (!check_init() || !check_init(v_shpm)) {return v_shpm;}
  
  if (transforms.unit.size() == 0) {
    print_warning("transform", "No transformations were defined.");
    return v_shpm;
  }
  
  if (index < 0 || index >= (int)transforms.unit.size()) {
    print_warning("transform", "Invalid index passed.");
    return v_shpm;
  }
  
  VecShpm result;
  for (uint i=0; i<(uint)v_shpm.size(); i++) {
  result.push_back(intern_transform(v_shpm[i], transforms.unit[index]));
  }
  
  return result;
  
}

Shpm Shc::create_weighting_hemi() {
  
  if (!check_init()) {return Shpm();}
  
  // calculate coefficients of weighting function (cos(theta)+1)/2  
  Coef coef(n_bands*n_bands);
  coef.setZero();
  coef(0) = sqrt(M_PI/1.0) / (2.0*M_PI);
  coef(2) = sqrt(M_PI/3.0) / (2.0*M_PI);
  
  Shpm result = coef2shpm(coef, HEMI_RMN, 2);
  
  return result;
  
}

Xyz Shc::compass(VecShpm& v_shpm1, VecShpm& v_shpm2) {
   
  VecShpm s1 = v_shpm1;
  VecShpm s2 = v_shpm2;
  VecShpm m1;
  VecShpm m2;

 return compass(s1, s2, m1, m2);
  
}
Xyz Shc::compass(VecShpm& v_shpm1, VecShpm& v_shpm2, VecShpm& v_mask1, VecShpm& v_mask2) {
  
  if (!check_init() || !check_init(v_shpm1)|| !check_init(v_shpm2) || !check_init(v_mask1) || !check_init(v_mask2)) {return Xyz();}
  
  if (v_mask1.size() > 0 || v_mask2.size() > 0) {
    if (equal(v_shpm1.size(), v_shpm2.size(), v_mask1.size(), v_mask2.size()) == false) {
      print_warning("compass", "Input vectors (shpm, masks) have different sizes");
      return Xyz();
    }
  } else {
    if (equal(v_shpm1.size(), v_shpm2.size()) == false) {
      print_warning("compass", "Input vectors have different sizes");
      return Xyz();
    }
  }
  
  if (v_mask1.size() > 0 && n_bands_CG == 0) {
    print_warning("compass", "To use weighting functions set the number of CG bands to a value greater than zero");
  }
  
  int n = v_shpm1.size();
  if (n != lut_translate.n_slices) {
    print_warning("compass", "Size of input vectors has to equal the number of slices.");
    return Xyz();
  }
 
  VecShpm s1 = v_shpm1;
  VecShpm s2 = v_shpm2;
  VecShpm m1 = v_mask1;
  VecShpm m2 = v_mask2;
  
  Xyz xyz = intern_compass(v_shpm1, v_shpm2, m1, m2);
  
  return xyz;
  
}

MatrixReal Shc::product(MatrixReal& matrix1, MatrixReal &matrix2) {
  
  MatrixReal result;
  
  if (matrix1.cols() != matrix2.cols() || matrix1.rows() != matrix2.rows()) {
    print_warning_static("product", "matrix dimensions do not match.");
    return result;
  }
  
  result = matrix1.cwiseProduct(matrix2);
  return result;
  
}
MatrixReal Shc::rotate(MatrixReal& matrix, Xyz xyz) {
  return intern_warp(matrix, Coor3d(0,0,0), xyz, VISUAL);
}
MatrixReal Shc::warp(MatrixReal& matrix, Coor3d translation, e_trans_type trans_type) { 
  return intern_warp(matrix, translation, Xyz(0,0,0), trans_type);
}

MatrixReal Shc::resize(MatrixReal& matrix, int width, int height) {
  
  int w = matrix.cols();
  int h = matrix.rows();
  
  MatrixReal result(height, width);
  
  if (h == 0 || w == 0) {
    print_warning_static("resize", "matrix is empty");
    result.setZero();
    return result;
  }
  if (width <= 0 || height <= 0) {
    print_warning_static("resize", "invalid size");
    result.setZero();
    return result;
  }
  
  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      float my = (float)(y)/(float)(height) * h;
      float mx = (float)(x)/(float)(width) * w;
      result(y,x) = filter_bilinear_mat(matrix, h, w, my, mx);
    }
  }
  
  return result;
  
}
MatrixReal Shc::convolution(MatrixReal& matrix, MatrixReal &kernel) {
  
  int mI = matrix.rows(); int nI = matrix.cols();
  int mK = kernel.rows(); int nK = kernel.cols();
  
  MatrixReal M = MatrixReal::Zero(mI,nI);
  
  if (mK%2 == 0 || nK%2 == 0) {
    print_warning_static("convolution", "Only kernel matrices with odd dimensions can be used");
    return M;
  }
  
  int mS = (mK-1)/2;
  int nS = (nK-1)/2;
  
  for (int row=mS; row<mI-mS; row++) {
    for (int col=nS; col<nI-nS; col++) {   
      M.coeffRef(row,col) = (matrix.block(row-mS,col-nS,mK,nK).cwiseProduct(kernel)).sum();
    }
  }

  for (int row=mS; row<mI-mS; row++) {
    for (int col=0; col<nS; col++) {
      float nR = nS-col;
      M.coeffRef(row,col)  = (matrix.block(row-mS,0,mK,nK-nR).cwiseProduct(kernel.block(0,nR,mK,nK-nR))).sum();
      M.coeffRef(row,col) += (matrix.block(row-mS,nI-nR,mK,nR).cwiseProduct(kernel.block(0,0,mK,nR))).sum();
    }
    for (int col=nI-nS; col<nI; col++) {
      float nR = col-nS;
      M.coeffRef(row,col)  = (matrix.block(row-mS,nR,mK,nI-nR).cwiseProduct(kernel.block(0,0,mK,nI-nR))).sum();
      M.coeffRef(row,col) += (matrix.block(row-mS,0,mK,nK-nI+nR).cwiseProduct(kernel.block(0,nI-nR,mK,nK-nI+nR))).sum();
    }
  }

  for (int col=nS; col<nI-nS; col++) {   
    for (int row=0; row<mS; row++) {
      float mR = mS-row;
      M.coeffRef(row,col)  = (matrix.block(0,col-nS,mK-mR,nK).cwiseProduct(kernel.block(mR,0,mK-mR,nK))).sum();
      M.coeffRef(row,col) += (matrix.block(mI-mR,col-nS,mR,nK).cwiseProduct(kernel.block(0,0,mR,nK))).sum();
    }
    for (int row=mI-mS; row<mI; row++) {
      float mR = row-mS;
      M.coeffRef(row,col)  = (matrix.block(mR,col-nS,mI-mR,nK).cwiseProduct(kernel.block(0,0,mI-mR,nK))).sum();
      M.coeffRef(row,col) += (matrix.block(0,col-nS,mK-mI+mR,nK).cwiseProduct(kernel.block(mI-mR,0,mK-mI+mR,nK))).sum();
    }
  }
  
  return M;
  
}
MatrixReal Shc::convolution(MatrixReal& matrix, Matrix3f& kernel) {
    
  int mI = matrix.rows(); int nI = matrix.cols();
  
  MatrixReal M = MatrixReal::Zero(mI,nI);
  
  for (int row=1; row<mI-1; row++) {
    for (int col=1; col<nI-1; col++) {   
      M.coeffRef(row,col) = (matrix.block<3,3>(row-1,col-1).cwiseProduct(kernel)).sum();
    }
  }
  
  for (int row=1; row<mI-1; row++) {
    M.coeffRef(row,0)  = (matrix.block<3,2>(row-1,0).cwiseProduct(kernel.block<3,2>(0,1))).sum();
    M.coeffRef(row,0) += (matrix.block<3,1>(row-1,nI-1).cwiseProduct(kernel.block<3,1>(0,0))).sum();
    M.coeffRef(row,nI-1)  = (matrix.block<3,2>(row-1,nI-2).cwiseProduct(kernel.block<3,2>(0,0))).sum();
    M.coeffRef(row,nI-1) += (matrix.block<3,1>(row-1,0).cwiseProduct(kernel.block<3,1>(0,2))).sum();
  }
  
  M.block(   0,0,1,nI) = matrix.block(   0,0,1,nI);
  M.block(mI-1,0,1,nI) = matrix.block(mI-1,0,1,nI);
  
  return M;
  
}
MatrixReal Shc::convolution(MatrixReal& matrix, e_convolution_type convolution_type) {
  
  Matrix3f kh(3,3);
  Matrix3f kv(3,3);
  
  bool hv = false;
  
  switch (convolution_type) {
    case SOBEL:
      kh << 1, 0, -1,
            2, 0, -2,
            1, 0, -1;
      kv = kh.transpose();
      hv = true;
      break;    
    case PREWITT:
      kh << -1, 0,  1,
            -1, 0,  1,
            -1, 0,  1;
      kv = kh.transpose();
      hv = true;
      break;
    case EDGE:
      kh <<  0,  1,  0,
             1, -4,  1,
             0,  1,  0;
      break;
    case LAPLACE:
      kh <<  1,  2,  1,
             2,-12,  2,
             1,  2,  1;
      break;
    case BLUR_BOX:
      kh = Matrix3f::Ones();
      kh *= 1.0/9.0;
      break;
    case BLUR_GAUSSIAN:
      kh <<  1,  2,  1,
             2,  4,  2,
             1,  2,  1;
      kh *= 1.0/16.0;
      break;
  }
  
  MatrixReal resh, resv, res;
  if (hv == true) {
    resh = convolution(matrix, kh);
    resv = convolution(matrix, kv);
    res  = resh.array().abs() + resv.array().abs();
  } else {
    res = (convolution(matrix, kh)).array().abs();
  }
  
  return res;
  
}
MatrixReal Shc::histogram_equalization(MatrixReal& matrix, int bins_x, int bins_y, bool panoramic) {
  
  MatrixReal result;
  
  if (bins_y <= 0 || bins_x <= 0) {
    print_warning_static("histogramEqualization","The bins have to be 1 (histogram equalization) or greater (local histogram equalization)");
    return result;
  }
  
  int m = matrix.rows();
  int n = matrix.cols();
  
  if (bins_y == 1 && bins_x == 1) {
      
    // histogram equalization
    VectorReal hist = intern_calcHistogram(matrix, 0, M_PI, true, panoramic);

    result.resize(m,n);  
    for (int i=0; i<m*n; i++) {
      result(i) = intern_histogramEqualization_calcEntry(hist, matrix(i));
    }
    
  } else {
    
    // local histogram equalization   
    float yi = (float)(m) / (bins_y);
    float xi = (float)(n) / (bins_x);
    
    // create histogram for each block
    vector<VectorReal> v_hist;
    v_hist.reserve(bins_y*bins_x);
    
    for (int y=0; y<bins_y; y++) {
      for (int x=0; x<bins_x; x++) {

        int y0 = y*yi; int y1 = (y+1)*yi;
        int x0 = x*xi; int x1 = (x+1)*xi;

        MatrixReal block = matrix.block(y0, x0, y1-y0, x1-x0);
        VectorReal hist = intern_calcHistogram(block, float(y0)/float(m)*M_PI, float(y1)/float(m)*M_PI, true, panoramic);

        v_hist.push_back(hist);
        
      }
    }
    
    result.resize(m,n);
    for (int y=0; y<m; y++) {
      for (int x=0; x<n; x++) {
        
        float val = matrix(y,x);
        int ix0, iy0, ix1, iy1;
        float py, px;
        
        if (x >= xi/2 && x < n-xi/2) {
          float xx = x-xi/2.0f;
          ix0 = xx / xi;
          ix1 = ix0+1;
          px = fmod(xx / xi, 1.0f) ;
        } else {
          ix0 = bins_x-1;
          ix1 = 0;
          if (x < xi/2) {
            px = 0.5+x/xi;
          } else {
            px = 0.5-(n-x)/xi;
          }
        }

        if (y >= yi/2 && y < m-yi/2) {
          float yy = y-yi/2.0f;
          iy0 = yy / yi;
          iy1 = iy0+1;
          py = fmod(yy / yi, 1.0f) ;
        } else {
           if (y < yi/2) {
             iy0 = 0;
             iy1 = 0;
             py = 1;
           } else {
             iy0 = bins_y-1;
             iy1 = bins_y-1;
             py = 1;
           }
        }
 
        float f00 = intern_histogramEqualization_calcEntry(v_hist[iy0*bins_x+ix0], val);
        float f01 = intern_histogramEqualization_calcEntry(v_hist[iy0*bins_x+ix1], val);
        float f10 = intern_histogramEqualization_calcEntry(v_hist[iy1*bins_x+ix0], val);
        float f11 = intern_histogramEqualization_calcEntry(v_hist[iy1*bins_x+ix1], val);
        
        result(y,x) = filter_bilinear_calc(f00, f10, f01, f11, py, px);

      }
    }
    
  }
  
  return result;
  
}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
MatrixReal Shc::blur(MatrixReal& matrix, Coor3d axis, float sigma, bool inverse) {
  
  float an = axis.norm();
  if (an == 0) {
    print_warning_static("blur", "axis has length zero (null vector)");
    return matrix;
  } else {
    axis = axis.array() / an;
  }
  
  if (sigma <= 0) {
    if (sigma < 0) {
      print_warning_static("blur", "sigma has to be greater equal zero");
    }
    return matrix;
  }
  
  float s = sigma;
  int w = matrix.cols();
  int h = matrix.rows();
  float c = (float)(w)/2.0f;
  
  if (w == 0 || h == 0) {
    print_warning_static("blur", "matrix is empty");
    return matrix;
  }
  
  if (w%2 == 1) {
    print_warning_static("blur", "only matrices with an even number of columns can be blurred");
    return matrix;
  }
  
  // align rotation axis with Z axis (blur is always applied around the Z axis)
  s_sphericalCoor sc = coor32sphericalCoor(Coor3d(axis));
    
  Coor3d raxis;
  raxis(0) = cos(sc.phi+M_PI/2.0f);
  raxis(1) = sin(sc.phi+M_PI/2.0f);
  raxis(2) = 0;
  float angle = -sc.theta;
  
  Axr axr(angle, raxis);
  Xyz xyz = axr2xyz(axr);
  
  MatrixReal mrot = rotate(matrix, xyz);
  MatrixReal result(h,w);
  
  // prepare kernel (blur)
  VectorXd kernel(w);
  for (int x=0; x<w; x++) {
    kernel(x) = 1.0f/sqrt(2.0f*M_PI*s*s)*exp(-((x-c)*(x-c))/(2.0f*s*s));
  }
    
  FFT<double> fft;
  VectorXcd fft_kernel = fft.fwd(kernel);

  VectorXcd fft_result;
  
  // calc FFT on each image row
  VectorReal row;
  for (int y=0; y<h; y++) {
    row = mrot.row(y);
    VectorXd row2 = row.cast<double>();
    if (y>0 || y<h-1) {

      VectorXcd fft_image = fft.fwd(row2);
      
      if (inverse == false) {
        fft_result = fft_image.array() * fft_kernel.array();
      } else {
        for (int i=0; i<w; i++) {
          std::complex<double> c = fft_kernel[i];
          double sup = 1e-3;
          if (c.imag() >= 0 && c.imag() < sup) {
            c.imag(sup);
          }
          if (c.imag() <= 0 && c.imag() > -sup) {
            c.imag(-sup);
          }
          if (c.real() >= 0 && c.real() < sup) {
            c.real(sup);
          }
          if (c.real() <= 0 && c.real() > -sup) {
            c.real(-sup);
          }
          fft_kernel[i] = c;
        }
        fft_result = fft_image.array() / fft_kernel.array();
      }
      
      row2 = fft.inv(fft_result);
      
    }
    result.row(y) = row2.real().cast<float>();
  }
  
  // reorder matrix (inverse FFT mixes indices up)
  MatrixReal ret(h,w);
  ret.block(0,0,h,w/2)     = result.block(0,w/2,h,w/2);
  ret.block(0,w/2,h,w-w/2) = result.block(0,0,h,w/2);
  
  // undo Z axis alignment
  ret = rotate(ret, transpose(xyz)); 
  
  return ret;
  
}
MatrixReal Shc::thresholding(MatrixReal& matrix, float threshold) {
  
  int n = matrix.cols() * matrix.rows();
  
  if (n == 0)
    return matrix;
  
  MatrixReal result(matrix.rows(), matrix.cols());  
  for (int i=0; i<n; i++) {
    if (matrix(i) < threshold) {
      result(i) = 0;
    } else {
      result(i) = 1;
    }
  }
  
  return result;
  
}
MatrixReal Shc::thresholding_otsu(MatrixReal& matrix, float lambda, bool panoramic) {
  
  int n = matrix.cols() * matrix.rows();
  
  if (n == 0)
    return matrix;
  
  VectorReal hist = intern_calcHistogram(matrix, 0, M_PI, false, panoramic);
  float t = intern_otsu(hist);
  int rt = round(t*255.0f);

  if (lambda > 0) {
    float sum   = 0;
    float count = 0;
    for (int i = rt; i<256; i++) {
      sum   += i * hist(i);
      count += hist(i);
    }
    float mean = sum/count;
    
    sum = 0;
    for (int i = rt; i<256; i++) {
      sum += hist[i] * (i-mean)*(i-mean);
    }
    float std = sqrt(1.0f / count * sum);

    t = (mean - lambda*std) / 255.0f;
  }
    
  MatrixReal result = thresholding(matrix, t);
  
  return result;
  
}
MatrixReal Shc::project_plane(MatrixReal& matrix, float angle, int size) {
  
  MatrixReal out(size, size);
  
  int w = matrix.cols();
  int h = matrix.rows();
  
  angle /= 2.0f;
  
  if (angle <= 0 || angle >= M_PI/2.0f) {
    print_warning_static("project_plane", "Invalid viewing angle (limited to 0-180deg).");
  }
  if (size <  1) {
    print_warning_static("project_plane", "Invalid size (must be greater zero).");
  }
  
  float dist = 0.5 / tan(angle);
  
  for (float x=0; x<size; x++) {
    for (float y=0; y<size; y++) {
      
      s_coor3 coor;
      coor.x = x/size - 0.5;
      coor.y = y/size - 0.5;
      coor.z = dist;
      
      s_sphericalCoor sc = coor32sphericalCoor(coor);
      s_coor2 coor2 = sphericalCoor2coor2(sc);
      
      out(y,x) = filter_bilinear_mat(matrix, h, w, coor2.y*matrix.rows(), coor2.x*matrix.cols());
      
    }
  }
  
  return out;
  
}
MatrixReal Shc::project_fisheye(MatrixReal& matrix, float angle, int size) {
  
  if (angle < 0 || angle > 2.0f*M_PI) {
    print_warning_static("project_fisheye", "invalid opening angle, choose a value between 0 and 2*PI");
    return matrix;
  }
  if (size <= 0) {
    print_warning_static("project_fisheye", "size has to be greater zero");
    return matrix;
  }
  
  MatrixReal out(size, size);
  
  float h = matrix.rows();
  float w = matrix.cols();
  
  float z = (size-1)/2;
  for (int y=0; y<size; y++) {
    for (int x=0; x<size; x++) {
      
      float yy = (y-z)/z;
      float xx = (x-z)/z;
      float rr = sqrt(xx*xx+yy*yy);
      
      if (rr > 1) {
        out(y,x) = 0;
        continue;
      }
      
      float p = -(atan2(yy, xx)) + M_PI;
      if (p < 0)
        p+=2*M_PI;
      float t = rr * angle/2.0f;
      
      float st = sin(t);
      float weight;
      
      if (st == 0) {
        weight = w;
      } else {
        weight = 1.0f / st;
      }
      if (weight > 5) {
        weight = 5;
      }
      
      float cy = t/M_PI*h;
      float cx = p/(2.0f*M_PI)*w;
      
      out(y,x) = filter_bilinear_mat(matrix, h, w, cy, cx);
      
    }
  }
  
  return out;
  
}
MatrixReal Shc::normalize(MatrixReal& matrix) {
   
  int w = matrix.cols();
  int h = matrix.rows();
  
  MatrixReal result(h,w);
  
  float minValue = matrix.minCoeff();
  float maxValue = matrix.maxCoeff();
  
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      result(y,x) = (matrix(y,x) - minValue) / (maxValue - minValue);
    }
  }
  
  return result;
  
}
Surf Shc::normalize(Surf& surf) {
  
  if (!check_init()) {return surf;}
  
  int n = surf.size();
  if (n != lut_surf.n_angles) {
    print_warning_static("normalize", "the size of surf does not match the number of sample points on the surface");
  }
  
  float minValue = surf.minCoeff();
  float maxValue = surf.maxCoeff();
  for (int i=0; i<n; i++) {
    surf(i) = (surf(i) - minValue) / (maxValue - minValue);
  }
  
  return surf;
  
}

