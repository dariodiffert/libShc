
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_tangentDistance:
// This example shows how the tangent distance can be used to increase the invariance against arbitrary transformations.
// In this example two panoramic images with a big translational difference are loaded and rotationally misaligned.
// The rotational misalignment is shown before any correction and afterwards with tangent distance (with translation matrices) on and off.
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
  // set up the rotations internally used
  shx.init_rotations_cone(  16.0*M_PI/180.0,  64.0*M_PI/180.0);
  shx.init_rotations_cone(   8.0*M_PI/180.0,  16.0*M_PI/180.0);
  shx.init_rotations_cone(   4.0*M_PI/180.0,   8.0*M_PI/180.0);
  shx.init_rotations_cone(   2.0*M_PI/180.0,   4.0*M_PI/180.0);
  shx.init_rotations_cone(   1.0*M_PI/180.0,   2.0*M_PI/180.0);
  // we need to initialize it
  shx.init();

  // Only print important things to the console
  shx.set_print_level(WARNINGS);
  
  // Set some rotational misalignment which should be applied to our panoramic images
  Xyz rot1(10, 0,10,true);
  Xyz rot2( 0,10, 0,true);

  // precalculating the transformation matrices for the tangent distance
  cout << "precalculating tangent distance matrices ... " << flush;
  shx.add_tangent_distance_translation(AXIS_X, 0.1);
  shx.add_tangent_distance_translation(AXIS_Y, 0.1);
  shx.add_tangent_distance_translation(AXIS_Z, 0.1);
  cout << "done!" << endl;
  
  // load panoramic images
  MatrixReal mcv, mss;
  mcv = shx.load_matrix("material/tangentDistance/cv.bmp");
  mss = shx.load_matrix("material/tangentDistance/ss.bmp");
  
  mcv = shx.rotate(mcv, rot1);
  mss = shx.rotate(mss, rot2);
  
  // Fourier transform the panoramic images
  Shpm cv = shx.matrix2shpm(mcv);
  Shpm ss = shx.matrix2shpm(mss);
  
  // print rotation difference before the compass is applied  
  cout << "Rotational difference before compass: " << shx.angular_difference(rot1, rot2)*180.0/M_PI << " degree" << endl;
  
  // print rotation difference after the compass is applied (tangent distance OFF)
  shx.set_tangent_distance(OFF);  
  Xyz xyz1 = shx.compass(cv, ss);
  cout << "Rotational difference after compass (tangent distance OFF): " << shx.compass_evaluate(rot1, rot2, xyz1)*180.0/M_PI << " degree" << endl;
    
  // print rotation difference after the compass is applied (tangent distance ON)
  // the dual sided tangent distance can not be used with the compass since it would slow down the computation too much
  shx.set_tangent_distance(ONESIDED);
  Xyz xyz2 = shx.compass(cv, ss);
  cout << "Rotational difference after compass (tangent distance ON): " << shx.compass_evaluate(rot1, rot2, xyz2)*180.0/M_PI << " degree" << endl;
  
  shx.show({mcv, mss}, "Current view - Snapshot");
  
  return 0;
  
}








