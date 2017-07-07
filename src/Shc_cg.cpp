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
 
Shc::s_cgCoef Shc::cg_coefficients(int l, int m, int l1, int l2) {
  
  s_cgCoef result;
  
  float sum = 0; 
  float m1 = (m-l1-l2+abs(l1-l2+m))/2;
  int   n  = (m+l1+l2-abs(l1-l2-m))/2-m1+1;
  
  VectorReal B(2*n+1);
  B.setZero();
  VectorReal C(n+1+1);
  C.setZero();
  
  for (int x=n-1; x>=1; x--) {
    B(2*x) = l1*(l1+1)+l2*(l2+1)+2*(m1+x)*(m-m1-x)-l*(l+1);
    C(n) = 1;
    B(2*x-1) = sqrt((l1*(l1+1)-(m1+x)*(m1+x-1))*(l2*(l2+1)-(m-m1-x)*(m-m1-x+1)));
    C(x) = -(B(2*x)*C(x+1)+B(2*x+1)*C(x+2))/B(2*x-1);
    sum = sum + C(x)*C(x);
  }
  
  C(n) = sqrt(1.0/(sum+1));

  for (int x=n-1; x>=1; x--) {
    C(x) = -(B(2*x)*C(x+1)+B(2*x+1)*C(x+2))/B(2*x-1);
    sum = sum + C(x)*C(x);
  }
  
  
  MatrixReal MM(n,2);
  MM.setZero();
  
  int count = 0;
  for (int m1=-l1; m1 <= l1; m1++) {
    if (abs(m-m1) <= l2) {
      MM(count,0) = m1;
      MM(count,1) = m-m1;
      count++;
    }
  }
      
  result.cgIndices.resize(n,6);
  result.cgCoef.resize(n);
  result.column.resize(n);
  result.row.resize(n);
  
  int row;
  int column;
  for (int i=0; i<n; i++) {
    result.cgIndices(i,0) = l;
    result.cgIndices(i,1) = m;
    result.cgIndices(i,2) = l1;
    result.cgIndices(i,3) = MM(i,0);
    result.cgIndices(i,4) = l2;
    result.cgIndices(i,5) = MM(i,1);
    result.cgCoef(i)      = C(i+1);
    cg_indexConversion(l,m,l1,MM(i,0),l2,MM(i,1),row,column);
    result.row(i)    = row;
    result.column(i) = column;
  }
  result.n = n;
  

  return result;
  
}
void Shc::cg_indexConversion(int l, int m, int l1, int m1, int l2, int m2, int &row, int &column) {

    row = (l1+m1) * (2*l2+1) + l2 + m2;
    
    column = l+m;
    if (l > l2-l1) {
        column += l*l - (l2-l1)*(l2-l1); // Indizierung ueber summe 2*i+1 von a bis l-1
    }
    
}
MatrixRealSparse Shc::cg_createMatrix(int l1, int l2) {
  
  int n = (2*l1+1) * (2*l2+1);
  MatrixRealSparse result(n,n);
  
  vector<TD> list;
  
  s_cgCoef cg;
  
  for (int l=abs(l2-l1); l<=l2+l1; l++) {
    for (int m=-l; m<=l; m++) {

      cg = cg_coefficients(l,m,l1,l2);
      
      for (int k=0; k<cg.n; k++) {
        
        if (abs(cg.cgCoef(k)) >= tolerances.cg) { // !!! cutting all entries which are "nearly" zero or zero
          list.push_back(TD(cg.row(k),cg.column(k),cg.cgCoef(k)));
        }
      }
      
    }
  }

  result.setFromTriplets(list.begin(),list.end());
  
  return result;
  
}
MatrixRealSparse Shc::cg_createMatrixReal(int l1, int l2) {
  // the resulting matrix is complex, however ignoring the complex part does not affect any results, they cancel out.
  
  int n = (2*l1+1) * (2*l2+1);
  MatrixComplexSparse result(n,n);

  MatrixComplexSparse kp(n,n);
  KroneckerProductSparse<MatrixComplexSparse,MatrixComplexSparse>(lut_tm[l1],lut_tm[l2]).evalTo(kp);

  MatrixComplexSparse ds(n,n);
  ds = directSum(lut_tm,abs(l2-l1),l2+l1);
  
  result = kp.transpose() * lut_cg[l1][l2].cast<complex<float>>() * ds.conjugate();
  
  return result.real();
  
}
MatrixRealSparse Shc::cg_createMatrixRealCoupling(int l1, int l2) {
  // normalized refers to the coupling matrix
  
  float n = (2*l1+1) * (2*l2+1);
  MatrixRealSparse result(n,n);
  MatrixRealSparse f(n,n);
  
  vector<TD> list;
  float f_n, f_0, ff;
  
  int count=0;
  for (int l=abs(l2-l1); l<=l2+l1; l++) {
    
    f_n = sqrt( n / (2*l+1) );
    
    int r = (n-1)/2;
    int c = l*l-(l2-l1)*(l2-l1)+l;
    f_0 = lut_cg[l1][l2].coeff(r,c);
    ff = f_n*f_0;
    
    for (int i=0; i<2*l+1; i++) {
      if (abs(ff) >= tolerances.cg) { // !!! cutting all entries which are "nearly" zero or zero
        list.push_back( TD(count, count, ff) );
      }
      count++;
    }
  }
  
  f.setFromTriplets(list.begin(),list.end());
  
  result = lut_cg_real[l1][l2] * f;
  
  return result;
  
}
MatrixComplexSparse Shc::createTransformationMatrix(int l) {
  
  MatrixComplexSparse result(2*l+1,2*l+1);
  
//   if (l >= n_bands) // create dummy matrices for non reachable indices i=k+l
//     return result;
  
  vector<TC> list; 

  list.push_back(TC(l,l,complex<float>(sqrt(2),0)));
  
  for (int m=1; m<=l; m++) {
    list.push_back(TC(l+m,l+m,complex<float>(pow(-1,m),0)));
    list.push_back(TC(l+m,l-m,complex<float>(1,0)));
    list.push_back(TC(l-m,l+m,complex<float>(0,pow(-1,m-1))));
    list.push_back(TC(l-m,l-m,complex<float>(0,1)));
  }
  
  result.setFromTriplets(list.begin(),list.end());
  result = result / sqrt(2);
  result = result.adjoint().eval();

  return result;
  
}
MatrixComplexSparse Shc::directSum(vector<MatrixComplexSparse> &in, int start, int stop) { 
  // creates the directsum of the vector of matrices, in the given range (including start and stop)
  
  int dr=0;
  int dc=0;
  for (int i=start; i<=stop; i++) {
    dr += in[i].rows();
    dc += in[i].cols();
  }
  
  MatrixComplexSparse result(dr,dc);
  
  vector<TC> list;
  
  int cr = 0;
  int cc = 0;
  int mA; int nA;
  for (int i=start; i<=stop; i++) {

    for (int j=0; j<in[i].outerSize(); j++) {
      for (MatrixComplexSparse::InnerIterator it(in[i],j); it; ++it) {
        
        mA = it.row();
        nA = it.col();
        
        list.push_back( TC(cr+mA,cc+nA,it.value()) );
        
      }
    }
    
    cr += in[i].rows();
    cc += in[i].cols();
  }
  
  result.setFromTriplets(list.begin(),list.end());
  
  return result;
}
