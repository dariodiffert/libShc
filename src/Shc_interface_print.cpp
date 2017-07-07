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

void Shc::print(Xyz xyz) {
  
  stringstream ss;
  ss << fixed << setprecision(2);
  ss << "print rotation (xyz in deg): " << setw(7) << xyz.x*180/M_PI << ", " << setw(7) << xyz.y*180/M_PI << ", " << setw(7) << xyz.z*180/M_PI << endl;
  cout << ss.str() << endl;
  
}
void Shc::print(MatrixRotation r) {
  
  stringstream ss;
  ss << fixed << setprecision(2);
  ss << "print rotation (r): " << setw(5) << r(0,0) << " " << setw(5) << r(0,1) << " " << setw(5) << r(0,2) << endl;
  ss << "                    " << setw(5) << r(1,0) << " " << setw(5) << r(1,1) << " " << setw(5) << r(1,2) << endl;
  ss << "                    " << setw(5) << r(2,0) << " " << setw(5) << r(2,1) << " " << setw(5) << r(2,2) << endl;
  cout << ss.str() << endl;
  
}
void Shc::print(Axr axr) {
  
  stringstream ss;
  ss << fixed << setprecision(2);
  ss << "print rotation (angle in deg; axis): " << setw(7) << axr.angle*180.0f/M_PI << "; ";
  ss << setw(5) << axr.axis(0) << " " << setw(5) << axr.axis(1) << " " << setw(5) << axr.axis(2) << endl;
  cout << ss.str() << endl;
  
}
void Shc::print(Shpm shpm) {
  
  stringstream ss;
  ss << "print shpm: contin = ";
  switch (shpm.contin) {
    case FULL:        ss << "FULL, "; break;
    case HEMI_M:      ss << "HEMI_M, "; break;
    case HEMI_MN:     ss << "HEMI_MN, "; break;
    case HEMI_RM:     ss << "HEMI_RM, "; break;
    case HEMI_RMN:    ss << "HEMI_RMN, "; break;
  }
  ss << "l_max = " << shpm.l_max << ", #coef = " << shpm.coef.size();
  cout << ss.str() << endl;
  
}
void Shc::print(Pointcloud pointcloud) {
  
  stringstream ss;
  ss << "print pointcloud (size " << pointcloud.size() << "): " << endl;
  for (uint i=0; i<pointcloud.size(); i++) {
    ss << "(" << pointcloud[i][0] << ", " << pointcloud[i][1] << ", " << pointcloud[i][2] << ")" << endl;
  }
  ss << endl;
  cout << ss.str() << endl;
  
}
