
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_simple:
// A simple example which shows how to initialize the class, load a panoramic image from file, and Fourier transform it.
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
  // initialization phase: set a custom surface.
  // the surface information is stored in the file bee_eye.txt
  shx.init_surface("material/bee_eye.txt");
  // increase the number of bands to see more details
  shx.init_bands(20);
  shx.init_output(360,180);
  // we need to initialize it
  shx.init();
  
  // load an image and Fourier transform it
  Shpm shpm = shx.load_shpm("material/cv.bmp");
  
  // show the result
  // the black blob in the center of the image is the occiput (backside of the insects head) where no ommatidia are
  shx.show(shpm);
  
  return 0;
  
}








