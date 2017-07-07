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

void Shc::print(string message) {
  
  if (print_level == ALL) {
    cout << message << flush;
  }
  
}
void Shc::print(string prefix, float result, string suffix) {
  
  if (print_level == ALL) {
    cout << prefix << result << suffix << flush;
  }
  
}
void Shc::print_progress_start() {
  
  progress_current = 0; 
  
  if (print_level == ALL) {
    clock_gettime(CLOCK_MONOTONIC, &progress_time);
    cout << setw(3) << 0;
    cout << "%" << flush;
  }
  
}
void Shc::print_progress_update(float progress) {
  
  if (progress_current == -1) {
    return;
  }
  
  progress_current = progress;
  
  if (print_level == ALL) {
  
    struct timespec progress_time_current;
    clock_gettime(CLOCK_MONOTONIC, &progress_time_current);
    
    float ms = (progress_time_current.tv_sec * 1000.0 + progress_time_current.tv_nsec / 1000000.0) - (progress_time.tv_sec * 1000.0 + progress_time.tv_nsec / 1000000.0);
    
    if (ms > 500) {
      
      clock_gettime(CLOCK_MONOTONIC, &progress_time);
      cout << "\b\b\b\b" << setw(3) << (int)(progress_current*100);
      cout << "%" << flush;
      
    }
    
  }
  
}
void Shc::print_progress_stop() {
  
  progress_current = -1;
  
  if (print_level == ALL) {
    cout << "\b\b\b\b" << flush;
  }
  
}
void Shc::print_warning(string function, string warning) {
  
  if (print_level == ALL || print_level == WARNINGS) {
  
    if (progress_current != -1) {
      cout << "\b\b\b\b" << "progress timer was interrupted" << endl;
      print_progress_stop();
    }
    
    cout << "[Warning] in [" << function <<  "]: " << warning << endl;
    
  }
  
}
void Shc::print_warning_static(string function, string warning) {
  
  cout << "[Warning] in [" << function <<  " (static)]: " << warning << endl;
  
}