
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_compass:
// Loads two different views (currentview cv and snapshot ss) and aligns them using the compass.
// During initialization phase multiple rotations are initialized (coarse-to-fine search).
// -------------------------------------------
// -------------------------------------------

// check if example is executed from the correct directory
void check_path() {
  
  DIR* dir = opendir("material");
  if (dir) {
    closedir(dir);
  } else {
    cout << "Since local paths are used, examples have to be executed from their local directory" << endl;
    exit(0);
  }
    
}

int main() {
  
  check_path();
  
  // create an instance of Shx
  Shx shx;
  // initialization phase: set rotation parameters (coarse-to-fine).
  // we set them to cover ALL rotations, speed can be increased by limiting the search space.
  // the final resolution is (if the local minimum is not skipped due to the coarse-to-fine approach) 1.0 degree.
  shx.init_rotations_sphere(16.0*M_PI/180.0,  64.0*M_PI/180.0);
  shx.init_rotations_cone(   8.0*M_PI/180.0,  16.0*M_PI/180.0);
  shx.init_rotations_cone(   4.0*M_PI/180.0,   8.0*M_PI/180.0);
  shx.init_rotations_cone(   2.0*M_PI/180.0,   4.0*M_PI/180.0);
  shx.init_rotations_cone(   1.0*M_PI/180.0,   2.0*M_PI/180.0);
  shx.init_bands(15);
  // we need to initialize it
  shx.init();

  // load a currentview and snapshot (as matrices since we want to show them later)
  MatrixReal mcv = shx.load_matrix("material/cv.bmp");
  MatrixReal mss = shx.load_matrix("material/ss.bmp");
  
  // convert the matrices into the basis of real spherical harmonics
  Shpm cv = shx.matrix2shpm(mcv);
  Shpm ss = shx.matrix2shpm(mss);
  
  // use the visual compass to determine the rotational offset between both views
  Xyz xyz = shx.compass(cv, ss); 
  
  // print results
  shx.print(xyz);
  
  // apply the rotation to the current view
  Shpm rr = shx.rotate(cv, xyz);
  
  // show the results (Fourier transformed)
  shx.show({cv, rr, ss}, "current view - compass estimate - snapshot [Fourier transformed]", false);
  
  // show the results (original images)
  MatrixReal mrr = shx.rotate(mcv, xyz);
  shx.show({mcv, mrr, mss}, "current view - compass estimate - snapshot", true);
  
  return 0;
  
}








