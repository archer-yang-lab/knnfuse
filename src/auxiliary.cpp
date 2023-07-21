#include "gsf.h"

MatrixXd normaltheta (const MatrixXd& y, const MatrixXd& sigma, const MatrixXd& graph, const MatrixXd& wMtx, const MatrixXd& Eta, const MatrixXd& U){
  int k,d,i,j;
  MatrixXd Theta =  MatrixXd::Zero(D,K);
  MatrixXd A(D,D),
  B(D,1),
  C(D,1);
  MatrixXd I = MatrixXd::Zero(D,D);
  
  VectorXd wMtxSums(K);
  
  for(k = 0; k < K; k++) {
    wMtxSums(k) = wMtx.col(k).sum();
  }
  
  for(d=0 ; d < D; d++){
    I(d,d)=1 ;
  }
  //for(k = 0; k < K; k++) {
  k=0;
    A = wMtxSums[k]*sigma+graph.col(k).sum()*I;
    B = MatrixXd::Zero(D,1);
    for (i=0; i < n; i++){
      B = B + wMtx(i,k)*y.row(i);
    }
    C = MatrixXd::Zero(D,1);
    for (j = 0 ; j < K; j++){
      if (graph(k,j)==1){
        Rcpp::Rcout << "k=" << k << "j=" << j <<"\n";
        C = C + Eta.col(k+K*j)-U.col(k+K*j);
      }
    }
    //Theta.col(k) = A.inverse() * (B + C) ;
    Theta.col(k) = C ;
  //} 
  return Theta;
}

double etamax(const Matrix<double, 1, Dynamic>& z, double lambda){
  double normZ = z.norm(), u;
  u = 1-(1.0/normZ) * lambda;
  if( u>0.5) {
    return u;
  } else {
    return 0.5;
  }
}

VectorXd softThresholding(const VectorXd& z, double lambda){
  double c = 1 - (lambda / z.norm());

  if (c > 0) return c * z;
  else return VectorXd::Zero(D);
}

VectorXd scadUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  double normZ = z.norm();

  if(normZ <= (u+1) * lambda){
    return softThresholding(z, u * lambda);

  } else if( (u + 1) * lambda <= normZ && normZ < a * lambda) {
    return ((a - 1)/(a - u - 1)) * softThresholding(z, (a * u * lambda)/(a - 1));

  } else {
    return z;
  }
}

VectorXd mcpUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  double normZ = z.norm();

  if(normZ <= a * lambda){
    return (a/(a - u)) * softThresholding(z, u * lambda);

  } else {
    return z;
  }
}

// Note: in this case, "a" is the Adaptive Lasso weight. 
VectorXd adaptiveLassoUpdate(double u, const Matrix<double, 1, Dynamic>& z, double lambda, double a){
  return softThresholding(z, u * lambda * a);
}

VectorXd scadLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a){
  double scad, normEta = eta.norm();

  if (normEta <= lambda) {
    scad = lambda;

  } else if (lambda < normEta && normEta < a * lambda) {
    scad = (a * lambda - normEta) / (a - 1);

  } else {
    scad = 0;
  }

  return softThresholding(z, u * scad);
}

VectorXd mcpLLAUpdate(double u, const Matrix<double, 1, Dynamic>& z, const Matrix<double, 1, Dynamic>& eta, double lambda, double a){
  double mcp, normEta = eta.norm();

  if (normEta <= lambda * a) {
    mcp = (lambda - (normEta / a));

  } else {
    mcp = 0;
  }

  return softThresholding(z, u * mcp);
}






