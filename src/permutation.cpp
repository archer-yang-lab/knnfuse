#include "gsf.h"

bool isOrdered(int i, int j) {
  return (i < j);
}

double thetaDist(const VectorXd& theta1, const VectorXd& theta2) {
   double acc = 0.0;

   for (int i = 0; i < D; i++) {
     acc += pow((double)(theta1(i, 0) - theta2(i, 0)), 2);
   }

   return sqrt(acc);
 }

// Generates the symmetrix matrix (with 0 diagonal) of pairwise distances between the columns of theta.
MatrixXd getDistanceMatrix(const MatrixXd& theta) {
  MatrixXd distances(K, K);

  for (int k = 0; k < K; k++) {
    for (int j = 0; j < K; j++) {
      distances(k, j) = thetaDist(theta.col(j), theta.col(k));
    }
  }

  return distances;
}

// Linear search function.
bool find(std::vector<int> sigma, int j) {
  for (unsigned int i = 0; i < sigma.size(); i++) {
     if (sigma[i] == j) {
       return true;
     }
   }

   return false;
 }

// Generates the permutation alpha.
void alpha(const MatrixXd& theta, std::vector<int>& perm) {
  MatrixXd distances = getDistanceMatrix(theta);

   std::vector<int> sigma, tau;
   double maxEntry, tMinEntry, sMinEntry, tSum, sSum;
   int i, j, k, sResult, tResult;

   int argmaxInd[2];

   // Find the thetas which are most distant.
   maxEntry = -1;
   for (i = 0; i < distances.rows(); i++) {
     for (j = 0; j < distances.cols(); j++) {
       if (maxEntry < distances(i, j)) {
         maxEntry = distances(i, j);
         argmaxInd[0] = i;
         argmaxInd[1] = j;
       }
     }
  }

   sigma.push_back(argmaxInd[0]);
   tau.push_back(argmaxInd[1]);

   // Inductively move towards the nearest neighbor.
   for (k = 1; k < K; k++) {
     tMinEntry = sMinEntry = INT_MAX;
     sResult = -1;
     tResult = -1;

     for (j = 0; j < K; j++) {
       if (!find(sigma, j) && distances(j, sigma[k - 1]) <= sMinEntry) {
         sMinEntry = distances(j, sigma[k - 1]);
         sResult = j;
       }

       if (!find(tau, j) && distances(j, tau[k - 1]) <= tMinEntry) {
         tMinEntry = distances(j, tau[k - 1]);
         tResult = j;
       }
     }

     sigma.push_back(sResult);
     tau.push_back(tResult);
   }

   // Determine whether tau or sigma defines the shortest path, and choose it
   // as the permutation alpha.
   tSum = sSum = 0;
   for (k = 0; k < K; k++) {
     tSum += distances(k, tau[k]);
     sSum += distances(k, sigma[k]);
   }

   if (tSum < sSum) {
     for (k = 0; k < K; k++) {
       perm[k] = tau[k];
     }

   }
   else {
     for (k = 0; k < K; k++) {
       perm[k] = sigma[k];
     }
  }
 }

 // Reorders the columns of theta with respect to the permutation alpha.
 MatrixXd reorderTheta(const MatrixXd& theta) {
   MatrixXd result(D, K);
   std::vector<int> perm(K);

   alpha(theta, perm);

   for (int k = 0; k < K; k++) {
     result.col(k) = theta.col(perm[k]);
   }

   return result;
 }

Psi reorderResult(const Psi& psi) {
  MatrixXd thetaResult(D, K);
  VectorXd piiResult(K);
  std::vector<int> perm(K);
  Psi newPsi;

  alpha(psi.theta, perm);

  for(int k = 0; k < K; k++){
    thetaResult.col(k) = psi.theta.col(perm[k]);
    piiResult(k)       = psi.pii(perm[k]);
  }

  newPsi.theta = thetaResult;
  newPsi.pii   = piiResult;
  newPsi.sigma = psi.sigma;

  return newPsi;
}

Rcpp::IntegerVector smallestKIndices(VectorXd vec, int K) {
  int n = vec.size();
  
  if (K <= 0 || K > n) {
    Rcpp::stop("Invalid value of K");
  }
  
  // Create a vector of indices from 1 to n
  Rcpp::IntegerVector indices = Rcpp::seq(1, n);
  
  // Sort the indices based on the corresponding vector elements
  std::sort(indices.begin(), indices.end(), [&vec](int a, int b) {
    return vec(a - 1) < vec(b - 1);
  });
  
  // Extract the first K indices
  Rcpp::IntegerVector result = Rcpp::head(indices, K);
  
  return result;
}


// Generates the graphs
//gsf
MatrixXd graphgsf(const MatrixXd& theta) {
  MatrixXd distances = getDistanceMatrix(theta);
  MatrixXd graph = MatrixXd::Zero(K,K);
  
  std::vector<int> sigma, tau;
  double maxEntry, tMinEntry, sMinEntry, tSum, sSum;
  int i, j, k, sResult, tResult;
  
  int argmaxInd[2];
  
  // Find the thetas which are most distant.
  maxEntry = -1;
  for (i = 0; i < distances.rows(); i++) {
    for (j = 0; j < distances.cols(); j++) {
      if (maxEntry < distances(i, j)) {
        maxEntry = distances(i, j);
        argmaxInd[0] = i;
        argmaxInd[1] = j;
      }
    }
  }
  
  sigma.push_back(argmaxInd[0]);
  tau.push_back(argmaxInd[1]);
  
  // Inductively move towards the nearest neighbor.
  for (k = 1; k < K; k++) {
    tMinEntry = sMinEntry = INT_MAX;
    sResult = -1;
    tResult = -1;
    
    for (j = 0; j < K; j++) {
      if (!find(sigma, j) && distances(j, sigma[k - 1]) <= sMinEntry) {
        sMinEntry = distances(j, sigma[k - 1]);
        sResult = j;
      }
      
      if (!find(tau, j) && distances(j, tau[k - 1]) <= tMinEntry) {
        tMinEntry = distances(j, tau[k - 1]);
        tResult = j;
      }
    }
    
    sigma.push_back(sResult);
    tau.push_back(tResult);
  }
  
  // Determine whether tau or sigma defines the shortest path, and choose it
  // as the permutation alpha.
  tSum = sSum = 0;
  /*   for (k = 0; k < K; k++) {
   tSum += distances(k, tau[k]);
   sSum += distances(k, sigma[k]);
  }
   */   
  for (k = 1; k < K; k++) {
    tSum += distances(tau[k-1], tau[k]);
    sSum += distances(sigma[k-1], sigma[k]);
  }
  
  if (tSum < sSum) {
    for (k = 1; k < K; k++) {
      graph(sigma[k-1],sigma[k]) = 1;
      graph(sigma[k],sigma[k-1]) = 1;
    }
    
  }
  else {
    for (k = 1; k < K; k++) {
      graph(tau[k-1],tau[k]) = 1;
      graph(tau[k],tau[k-1]) = 1;
    }
  }
  return graph;
}

/*1nn*/
MatrixXd graph1nn(const MatrixXd& theta) {
  int K = theta.cols();
  MatrixXd distances = getDistanceMatrix(theta);
  MatrixXd graph = MatrixXd::Zero(K,K);
  
  double sMinEntry;
  int j, k, sResult;
  
  // Inductively move towards the nearest neighbor.
  for (k = 0; k < K; k++) {
    sMinEntry = INT_MAX;
    sResult = -1;
    
    for(j=0; j < K; j++){
      if (!(j==k) && distances(k, j) <= sMinEntry) {
        sMinEntry = distances(k,j);
        sResult = j;
      }
    }
    graph(k,sResult) = 1;
    graph(sResult,k) = 1;
    //Rcpp::Rcout << "k-th row" << graph.row(k) << "\n";
  }
  
  return graph;
}

/*mnn*/
MatrixXd graphmnn(const MatrixXd& theta, int m) {
  int K = theta.cols();
  MatrixXd distances = getDistanceMatrix(theta);
  MatrixXd graph = MatrixXd::Zero(K,K);
  Rcpp::IntegerVector nearestm(m, 0);
  
  int j, k;
  
  // Inductively move towards the nearest neighbor.
  for (k = 0; k < K; k++) {
    nearestm = smallestKIndices(distances.col(k),m+1);
    for (j=1; j < m+1; j++){
      graph(k,(nearestm(j)-1)) = 1;
      graph((nearestm(j)-1),k) = 1;
    }
  }
  return graph;
}


//naive
MatrixXd graphnaive(const MatrixXd& theta) {
  MatrixXd graph = MatrixXd::Ones(K,K);
  for (int j = 0; j < K; j++){
    graph(j,j) = 0;
  }
  return graph;
}

bool linSearch(const MatrixXd& target, const MatrixXd& list) {
  for (unsigned int i = 0; i < list.cols(); i++) {
    if (thetaDist(list.col(i), target) < 1e-6) {
      return true;
    }
   }

  return false;
}

int frequency(const MatrixXd& theta) {
	MatrixXd uniqueThetas;
	MatrixXd currentTheta(D, 1);

	for (int k = 0; k < K; k++) {
		currentTheta = theta.col(k);

      if (uniqueThetas.size() == 0 || !linSearch(currentTheta, uniqueThetas)) {
			uniqueThetas.conservativeResize(D, uniqueThetas.cols() + 1);
			uniqueThetas.block(0, uniqueThetas.cols() - 1, D, 1) = currentTheta;
		}
	}

	return uniqueThetas.cols();
}
