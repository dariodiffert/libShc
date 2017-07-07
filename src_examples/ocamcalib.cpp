
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_ocamcalib:
// loads a calib_results.txt file to unwrap an image using OCamCalib
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
  shx.init_bands(20);
  // load the calib_results.txt file in which the calibration parameters are stored to unwrap raw images
  shx.init_ocamcalib("material/ocamcalib/calib_results.txt", false, true, true);
  // we need to initialize it
  shx.init();
  
  MatrixReal m_in = shx.load_matrix("material/ocamcalib/panoramic.bmp");
  MatrixReal m_ocam = shx.ocamcalib2matrix(m_in);
  
  shx.show(m_in, "input image", false);
  shx.show(m_ocam, "panoramic image", true);
  
  return 0;
  
}








