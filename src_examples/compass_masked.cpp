
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_compass_masked:
// We enhance the example_compass by adding a mask to deal with corrupted images.
// In this example we assume that a part of the images is invalid (e.g. the camera setup is visible) and mask it out (see loaded the input images).
// It can be seen that, without using masks, the invalid area corrupts the result,
// while using masks can avoid this behaviour.
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
  // we want to use masks (which need products of spherical harmonics)
  // therefore we have to set the number of Clebsch-Gordan coefficients greater zero (depending on the accuracy)
  // since the mask have a simple form, a low number of bands should be sufficient
  shx.init_bands(20, 10);
  // we need to initialize it
  shx.init();

  // load two views, ss was recorded at a different location than ss and tilted.
  // additionally both images are corrupted
  MatrixReal mcv = shx.load_matrix("material/compass_masked/cv.bmp");
  MatrixReal mss = shx.load_matrix("material/compass_masked/ss.bmp");
  // now we load the mask (in this example the same mask for both images)
  MatrixReal mmm = shx.load_matrix("material/compass_masked/mask.bmp");
  
  // convert the views into the basis of real spherical harmonics
  Shpm cv = shx.matrix2shpm(mcv);
  Shpm ss = shx.matrix2shpm(mss);
  // since the mask has a simple structure, we use only a small number of bands to approximate it
  Shpm mm = shx.matrix2shpm(mmm, FULL, 5);
  
  // determine the rotation to align the views without masks
  Xyz xyz_simple = shx.compass(cv, ss);
  // determine the rotation to align the views using masks (both masks are the same, different masks for cv and ss can be passed!)
  Xyz xyz_masked = shx.compass(cv, ss, mm, mm);
  
  // print results
  cout << "best rotation found without masks: " << endl;
  shx.print(xyz_simple);
  cout << "best rotation found using masks: " << endl;
  shx.print(xyz_masked);
  
  // apply the rotation to cv
  Shpm rs = shx.rotate(cv, xyz_simple);
  Shpm rm = shx.rotate(cv, xyz_masked);
  
  // show the results (Fourier transformed)
  shx.show({cv, rs, rm, ss}, "Current view - Compass estimate (no mask) - Compass estimate (with mask) - Snapshot [Fourier transformed]", false);
  
  // show the results (original images)
  MatrixReal mrs = shx.rotate(mcv, xyz_simple);
  MatrixReal mrm = shx.rotate(mcv, xyz_masked);
  shx.show({mcv, mrs, mrm, mss}, "Current view - Compass estimate (no mask) - Compass estimate (with mask) - Snapshot", true);
  
  return 0;
  
}








