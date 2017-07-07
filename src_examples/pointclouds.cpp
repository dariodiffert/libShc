
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_pointclouds:
// here we load a 3D model (source: http://visionair.ge.imati.cnr.it/ontologies/shapes/viewgroup.jsp?id=670-fertility_-_watertight) and load it into a pointcloud.
// a second rotated version of the pointcloud is created and we want to find the rotation between both of them using the compass.
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
  // set the number of bands used
  shx.init_bands(15);
  // set up coarse-to-fine rotations for the compass
  shx.init_rotations_sphere(16.0*M_PI/180.0,  64.0*M_PI/180.0);
  shx.init_rotations_cone(   8.0*M_PI/180.0,  16.0*M_PI/180.0);
  shx.init_rotations_cone(   4.0*M_PI/180.0,   8.0*M_PI/180.0);
  shx.init_rotations_cone(   2.0*M_PI/180.0,   4.0*M_PI/180.0);
  shx.init_rotations_cone(   1.0*M_PI/180.0,   2.0*M_PI/180.0);
  // we can set the number of slices greater 1 to obtain more depth information of the structure
  shx.init_slices(1);
  shx.init_output(360,180);
  // we need to initialize it
  shx.init();
  

  // we load the pointcloud of a 3D scan of a sculpture 
  // Vertices: 143, source: http://visionair.ge.imati.cnr.it/ontologies/shapes/viewgroup.jsp?id=670-fertility_-_watertight
  Pointcloud pointcloud_cv = shx.load_pointcloud_auto("material/pointcloud.txt");
  
  // We create a second, rotated (in degree), version of the sculpture ...
  Xyz rot(10, 20, 30, true);
  Pointcloud pointcloud_ss = shx.rotate(pointcloud_cv, rot);
  
  // ... and Fourier transform them
  VecShpm vshpm_cv = shx.pointcloud2shpm(pointcloud_cv);
  VecShpm vshpm_ss = shx.pointcloud2shpm(pointcloud_ss);

  // determine rotation betwen both
  Xyz xyz = shx.compass(vshpm_cv, vshpm_ss);
  
  VecShpm vshpm_rr = shx.rotate(vshpm_cv, xyz);
  
  cout << "rotation applied to the pointcloud" << endl;
  shx.print(rot);
  cout << "rotation determined by the visual compass" << endl;
  shx.print(xyz);

  return 0;
  
}








