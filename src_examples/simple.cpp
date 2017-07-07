
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
  
  // set the number of bands used to compute the Fourier coefficient vector
  const int bands = 6;
  
  // create an instance of Shx
  Shx shx;
  // initialization phase: set the number of bands used
  shx.init_bands(bands);
  // we need to initialize it
  shx.init();
 
  // load an image and Fourier transform it. The data is stored in shpm
  Shpm shpm = shx.load_shpm("material/cv.bmp");
  
  // show some general information
  shx.print(shpm);
  
  // get the Fourier coefficient vector
  Coef coef = shx.shpm2coef(shpm);
  
  // just for convinience, get maximal number of bands which equals the number of bands = 6 used above
  int l_max = shpm.get_max_band();
  
  // print the Fourier coefficients
  cout << "The Fourier coefficients are: " << endl;
  int c=0;
  for (int l=0; l<l_max; l++) { 
    for (int m=-l; m<=l; m++) { 
      cout << "l = " << l << ", m = " << m << ", val = " << coef(c++) << endl;
    }
  }  
    
  return 0;
  
}








