
// include the Shc library (with Shx extension)
#include "Shx.h"
#include <dirent.h>

// setup namespaces
using namespace std; 
using namespace Eigen;
using namespace shc;

// -------------------------------------------
// -------------------------------------------
// example_seqSlam:
// The seqSlam algorithm has been implemented to identify the position of an agent (during the second transit) which drives two times the same track.
// The agent was mounted with a hemispherical camera pointing upwards, the images have been segmentated into sky and ground.
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

// SeqSLAM algorithm, adapted from Subokita (https://github.com/subokita/OpenSeqSLAM)
pair<int, float> seqSLAM(MatrixXf& M, int index, int matching_dist, float min_velocity, float max_velocity) {
  
  const float MAXVAL = 1e12;
  
  if (index < matching_dist / 2 + 1)
    return pair<int, float>(-1,0);
  
  if (index > M.cols() - matching_dist / 2 - 1)
    return pair<int, float>(-1,0);
  
  int move_min = round(min_velocity * matching_dist);
  int move_max = round(max_velocity * matching_dist);
  int y_max = M.rows();
  int n = move_max - move_min + 1;
  int m = matching_dist + 1;
  
  // create linspace for possible velocities
  VectorXf velocity(n);
  for(int i=0; i<n; i++) {
    velocity[i] = (move_min + i * 1.0) / matching_dist;
  }
  
  // Create incremental indices based on the previously calculated velocity
  MatrixXf increment_indices(n, m);
  for (int y=0; y<increment_indices.rows(); y++) {
    for (int x=0; x<increment_indices.cols(); x++) {
      increment_indices(y,x) = round(x * velocity(y));
    }
  }

  // Start trajectory
  int n_start = index - (matching_dist / 2);
  MatrixXf traj_start(n, m);
  for (int j=0; j<m; j++) {
    for (int i=0; i<n; i++) {
      traj_start(i,j) = (n_start + j - 1) * y_max;
    }
  }

  VectorXf score(y_max);
  
  // Perform the trajectory search to collect the scores
  for (int s=0; s<y_max; s++) {
    MatrixXf traj_dir = increment_indices + s * MatrixXf::Ones(n,m);
    
    for (int j=0; j<m; j++) {
      for (int i=0; i<n; i++) {
        if (traj_dir(i,j) > y_max) {
          traj_dir(i,j) = y_max;
        }
      }
    }
    
    MatrixXf traj_dest = traj_start + traj_dir;
    
    float min_sum = MAXVAL;
    for (int i=0; i<traj_dest.rows(); i++) {
        float sum = 0.0;
        for (int j=0; j<traj_dest.cols(); j++) {
          int dest = traj_dest(i, j);
          sum += M(dest % y_max, dest / y_max );
        }
        min_sum = min(min_sum, sum );
    }
    score(s) = min_sum;
  }
  
  // Find the lowest score
  int min_index = 0;
  float min_value = MAXVAL;
  for (int i=0; i<y_max; i++) {
    if (score(i) < min_value) {
      min_value = score(i);
      min_index = i;
    }
  }
  
  // ... now discard the region from where we found the lowest score
  for(int i = max(0, min_index - matching_dist / 2); i < min(y_max, min_index + matching_dist / 2); i++) {
    score[i] = MAXVAL;
  }
  
  // ... in order to find the second lowest score 
  float min_value_2 = MAXVAL;
  for (int i=0; i<y_max; i++) {
    if (score(i) < min_value_2) {
      min_value_2 = score(i);
    }
  }
  
  return pair<int, float> ( min_index + matching_dist / 2, min_value / min_value_2 );
    
}

// loads all images in the path (fullfilling the given name convention) and calculates the amplitude spectra in the basis of spherical harmonics
vector<VectorXf> image2AS(Shx &shx, string path) {
  
  vector<VectorXf> result;
  stringstream ss;
  Shpm shpm;
  
  int count = 1;
  while (true) {
    
    // read image
    stringstream number;
    number << std::setw(4) << std::setfill('0') << count;
    ss << path << "image-" << number.str() << ".png";
    shpm = shx.load_shpm(ss.str());
    ss.str("");
    
    // if the image could be read ...
    if (shpm.initialized() == false) {
      break;
    } else {
      
      // ... calculate the amplitude spectrum and push it on the result vector
      VectorXf amplitudeSpectrum = shx.get_feature(shpm, AS);
      result.push_back(amplitudeSpectrum);      
    }
    
    count++;
    
  }
  
  return result;
  
}

// calculate the difference matrix needed to do seqSLAM by comparing the amplitude spectra of all images of run1 with the images of run2
MatrixXf compute_diff_matrix(string run1, string run2) {
  
  // initialize shx (spherical harmonic stuff) instance
  Shx shx;
  shx.init_output(90, 45);
  shx.init_bands(40);
  shx.init_surface(4000);
  shx.init();

  // run through all images, transform them into the space of spherical harmonics, and calculate the amplitude spectrum
  vector<VectorXf> f1 = image2AS(shx, run1);
  vector<VectorXf> f2 = image2AS(shx, run2);

  // get number of loaded images
  int n1 = f1.size();
  int n2 = f2.size();
  
  // create difference matrix
  MatrixXf M(n2,n1);
  for (int i1=0; i1<n1; i1++) {
    for (int i2=0; i2<n2; i2++) {
      VectorXf diff = f1[i1] - f2[i2];
      M(i2,i1) = diff.norm();
    }
  }

  // transpose M to switch the matching direction: run1 -> run2 (transposed) or run2 -> run1 (not transposed).
  M.transposeInPlace();
  
  return M;
  
}

// start
int main() {
  
  check_path();
  
  // Path to folders containing preprocessed data (unfolded, thresholded)
  string run1 = "material/seqSlam/run1/";
  string run2 = "material/seqSlam/run2/";
  
  // Compute amplitude spectra for each image in each run.
  // Then create difference matrix, which stores the difference between any pair of images: diff(imageFromRun1, imageFromRun2)
  MatrixXf M;
  M = compute_diff_matrix(run1, run2);
  
  // find for each image i of run1 the best matching image result.first of run2. The confidence value is stored in result.second.
  for (int i=0; i<M.cols(); i++) {
    pair<int, float> result = seqSLAM(M, i, 10, 0.8, 1.2);
    cout << i << " " << result.first << " " << result.second << endl;
  }
  cout << "The output above shows which image of the first run (training) is matched with which image in the second run (test). The third value gives an estimate for the match quality." << endl;

}

