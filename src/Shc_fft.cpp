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

void Shc::Wfft::init(int n, int inverse) {

  this->n = n;
  this->inverse = inverse;
  
  void* mem_temp = NULL;
  size_t len_temp = 0;
  kiss_fftr_alloc(n, inverse, mem_temp, &len_temp);
  
  buffer.resize(len_temp);
  kiss_fftr_alloc(n, inverse, &buffer[0], &len_temp);
  
}  
Shc::Wfft::Wfft(int n, int inverse) {
  
  this->init(n, inverse);  
  
}  
Shc::Wfft::Wfft() {

  this->n = 0;
  this->inverse = 0;
  
}  
Shc::Wfft& Shc::Wfft::operator=(const Wfft& src) {
  
  if(&src == this) {
    return *this;
  }
      
  this->n = src.n;
  this->inverse = src.inverse;
  
  this->init(n, inverse);
  
  return *this;
  
}
Shc::Wfft::Wfft(const Wfft& src) {
  
  n = src.n;
  inverse = src.inverse;
  
  init(n, inverse);
  
}
  
void Shc::Wfft::compute(VectorReal& in, int offset, vector<kiss_fft_cpx>& out) {
  
  kiss_fftr((kiss_fftr_cfg)(&(buffer[0])), &(in[offset]), &(out[0]));  
  
}
void Shc::Wfft::compute(vector<kiss_fft_cpx>& in, VectorReal& out, int offset) {
  
  kiss_fftri((kiss_fftr_cfg)(&buffer[0]), &in[0], &out[offset]);
  
}

