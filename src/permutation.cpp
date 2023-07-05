#include "gsf.h"

bool isOrdered(int i, int j) {
  return (i < j);
}

double thetaDist(const VectorXd& theta1, const VectorXd& theta2) {
  int D = theta1.rows();  // Dimension of the parameter space.
  
   double acc = 0.0;

   for (int i = 0; i < D; i++) {
     acc += pow((double)(theta1(i, 0) - theta2(i, 0)), 2);
   }

   return sqrt(acc);
 }

// Generates the symmetrix matrix (with 0 diagonal) of pairwise distances between the columns of theta.
MatrixXd getDistanceMatrix(const MatrixXd& theta) {
  int K = theta.cols();  // Number of mixture components. 
  
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
  int K = theta.cols();  // Number of mixture components. 
  
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
   int K = theta.cols();  // Number of mixture components. 
   int D = theta.rows();  // Dimension of the parameter space.

   MatrixXd result(D, K);
   std::vector<int> perm(K);

   alpha(theta, perm);

   for (int k = 0; k < K; k++) {
     result.col(k) = theta.col(perm[k]);
   }

   return result;
 }

Psi reorderResult(const Psi& psi) {
  int K = psi.theta.cols();  // Number of mixture components. 
  int D = psi.theta.rows();  // Dimension of the parameter space.
  
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
// Generates the graphs
//gsf
MatrixXd graphgsf(const MatrixXd& theta) {
  MatrixXd distances = getDistanceMatrix(theta);
  int K = theta.cols();
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
  MatrixXd distances = getDistanceMatrix(theta);
  K = theta.cols();  // Number of mixture components. 
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
  }
  return graph;
}

//naive
MatrixXd graphnaive(const MatrixXd& theta) {
  K = theta.cols();  // Number of mixture components. 
  MatrixXd graph = MatrixXd::Ones(K,K);
  for (int j = 0; j < K; j++){
    graph(j,j) = 0;
  }
  return graph;
}

bool linSearch(const MatrixXd& target, const MatrixXd& list) {
  for (unsigned int i = 0; i < list.cols(); i++) {
    if (thetaDist(list.col(i), target) == 0) {
      return true;
    }
   }

  return false;
}

int frequency(const MatrixXd& theta) {
  int K = theta.cols();  // Number of mixture components. 
  int D = theta.rows();  // Dimension of the parameter space.
  
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
