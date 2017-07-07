
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_hemispherical_continuation:
// The upper hemisphere of a panoramic image can be mirrored to the lower hemisphere to fill in missing information.
// Using hemispherical continuation, the computation of most functions (e.g. the compass) is increased (commonly halved).
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
  shx.init_bands(30);
  shx.init_output(360,180);
  shx.init();

  // loads image with different hemispherical continuation RM
  Shpm shpm_full    = shx.load_shpm("material/cv.bmp");
  Shpm shpm_hemi_RM = shx.load_shpm("material/cv.bmp", HEMI_RM);
  
  // show the Fourier transformed image
  shx.show({shpm_full, shpm_hemi_RM}, "Full Fourier transform - Hemispherical Continuation");
  
  return 0;
  
}








