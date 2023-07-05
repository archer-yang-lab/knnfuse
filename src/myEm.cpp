#include "gsf.h"

int K, N, D, n, M, graphtype, maxadmm, maxPgd, maxRep, lambdaIter, modelIndex, penalty;
double epsilon, u, ck, tau1, a, delta,  
lambdaScale, uBound, H, alStart;

bool arbSigma, verbose;

MatrixXd upTransf, invEtaTransf, ones, transformedData;
std::vector<MatrixXd> yOuterProd; 
VectorXd lambdaVals, adaptiveLassoWeights;

MatrixXd (*gradB)(const MatrixXd&, const MatrixXd&);
double   (*b)(const VectorXd&, const MatrixXd&);
MatrixXd (*T)(const MatrixXd&);
double   (*density)(const Matrix<double, 1, Dynamic>&, const Matrix<double, Dynamic, 1>&, const MatrixXd&);
MatrixXd (*invTransf)(const MatrixXd&, const MatrixXd&);
MatrixXd (*Graphmat)(const MatrixXd&);
MatrixXd (*updateTheta)(const MatrixXd& , const MatrixXd&, const MatrixXd&, const MatrixXd&, const MatrixXd&, const MatrixXd& );
VectorXd (*updateEta)(double, const Matrix<double, 1, Dynamic>&, double, double);
VectorXd (*updateEtaLLA)(double, const Matrix<double, 1, Dynamic>&, const Matrix<double, 1, Dynamic>&, double, double);
bool     (*constrCheck)(const MatrixXd&);

Rcpp::List estimateSequence(const MatrixXd& y, const Psi& startingVals, const VectorXd& lambdaList);
Rcpp::List rbic(const MatrixXd&, const Psi&, const VectorXd&);
Psi mem(const MatrixXd&, const Psi&, double);
Matrix<double, Dynamic, Dynamic> wMatrix(const MatrixXd& y, const Psi&);
Psi mStep(const MatrixXd& , const Psi& , const MatrixXd& ,  const MatrixXd& , double );
MatrixXd admm(const MatrixXd& , const Psi&, const MatrixXd&, const MatrixXd&, double );
MatrixXd pgd(const MatrixXd&, const Psi&, const MatrixXd&, double lambda) ;
bool uCheck(double u, const MatrixXd& y, const MatrixXd& oldEta, const MatrixXd& newEta, const MatrixXd& sigma, const VectorXd& wMtxSums);


// [[Rcpp::export(.myEm)]]
extern "C" Rcpp::List myEm(SEXP argY, SEXP argGraphtype, SEXP argTheta, SEXP argSigma, 
                          SEXP argPii, SEXP argArbSigma, SEXP argM, 
                          SEXP argIndex, SEXP argMaxAdmm, SEXP argCk, SEXP argA, 
                          SEXP argPenalty, SEXP argLambdaVals, SEXP argEpsilon, 
                          SEXP argMaxRep, SEXP argDelta,
                          SEXP argVerbose){
  MatrixXd theta  = Rcpp::as<MatrixXd> (argTheta);      // Size D x K.
  //MatrixXd Graph  = Rcpp::as<MatrixXd> (argGraph);      // Size K x K.
  MatrixXd y      = Rcpp::as<MatrixXd> (argY);          // Size n x N.
  MatrixXd sigma  = Rcpp::as<MatrixXd> (argSigma);      // Size N x N.
  VectorXd pii    = Rcpp::as<VectorXd> (argPii);        // Size 1 x K.
  
  graphtype       = Rcpp::as<int> (argGraphtype);       //the graph to use in penalty
  arbSigma        = Rcpp::as<bool> (argArbSigma);       // Common unknown structure parameter check.
  ck              = Rcpp::as<double> (argCk);           // Parameter for penalty on the pi_k.
  epsilon         = Rcpp::as<double> (argEpsilon);      // Convergence criterion for EM algorithm.
  maxadmm         = Rcpp::as<int> (argMaxAdmm);
  //  maxPgd          = Rcpp::as<int> (argMaxPgd);          // Maximum number of repetetitions of the PGD algorithm.
  maxRep          = Rcpp::as<int> (argMaxRep);          // Maximum EM iterations.
  a               = Rcpp::as<double> (argA);            // Parameter for the SCAD or MCP penalty.  
  delta           = Rcpp::as<double> (argDelta);        // Convergence criterion for PGD algorithm.
  lambdaVals      = Rcpp::as<VectorXd> (argLambdaVals); // Values of lambda to be used.
  //  uBound          = Rcpp::as<double> (argUBound);       // Upper bound on the parameter u for the PGD algorithm. 
  verbose         = Rcpp::as<bool> (argVerbose);        // 
  M               = Rcpp::as<int> (argM);               // Number of trials for multinomial mixtures. 
  modelIndex      = Rcpp::as<int> (argIndex);           // Model index -- see below.
  //  H               = Rcpp::as<double> (argH);            // Hartigan lower bound for location-scale normal mixtures.
  penalty         = Rcpp::as<int> (argPenalty);                                   // H = 0 for other models. 
  
  // Set important constants.
  K = theta.cols();  // Number of mixture components. 
  D = theta.rows();  // Dimension of the parameter space.
  N = y.cols();      // Dimension of the sample space. 
  n = y.rows();      // Sample size.
  //set graph type 
  switch(graphtype){
  case 1: Graphmat = &graphnaive;    //naive pairwise 
  case 2: Graphmat = &graph1nn; //1nn
  case 3: Graphmat = &graphgsf; //gsf
    
  }
  
  
  // Initialize the set of parameters.
  Psi psi;
  psi.pii   = pii;
  psi.sigma = sigma;
  
  alStart = 0;
  
  switch(modelIndex) {
  
  // User selected a multivariate normal mixture in location. 
  // The covariance matrix is assumed to be common, but if hasFixedMatrix is false, 
  // it is unknown and estimated by the EM algorithm.
  case 1:   updateTheta = &normaltheta;
    gradB       = &gradBNormalLoc;
    b           = &bNormalLoc;
    T           = &tNormalLoc;
    density     = &densityNormalLoc;
    invTransf   = &invTransfNormalLoc;
    constrCheck = &constrCheckNormalLoc;            
    
    psi.theta   = transfNormalLoc(theta, sigma);
    
    for (int i = 0; i < n; i++) {
      yOuterProd.push_back((y.row(i).transpose())*(y.row(i)));
    }
    
    break;
    // User selected a multinomial mixture.
  case 3: gradB       = &gradBMultinomial;
    b           = &bMultinomial; 
    T           = &tMultinomial;
    density     = &densityMultinomial;
    invTransf   = &invTransfMultinomial;
    constrCheck = &constrCheckMultinomial;
    
    psi.theta   = transfMultinomial(theta);
    
    break;
    
    // User selected a Poisson mixture.
  case 5: gradB       = &gradBPoisson; 
    b           = &bPoisson;
    T           = &tPoisson;
    density     = &densityPoisson;
    invTransf   = &invTransfPoisson;
    constrCheck = &constrCheckPoisson;
    
    psi.theta   = transfPoisson(theta);
    
    break;
    
    // User selected a mixture of exponential distributions.
  case 6: gradB       = &gradBExponential; 
    b           = &bExponential;
    T           = &tExponential;
    density     = &densityExponential;
    invTransf   = &invTransfExponential;
    constrCheck = &constrCheckExponential;
    
    psi.theta   = transfExponential(theta);
    
    break;
  }
  
  transformedData = (*T)(y).transpose();
  Rcpp::List out = estimateSequence(y,psi, lambdaVals);
  
  return out; 
}

// The Modified EM (MEM) Algorithm wrapper.
Psi mem(const MatrixXd& y, const Psi& psi, double lambda) {
  MatrixXd wMtx, th, graph;
  Psi oldEstimate, newEstimate = psi;
  int counter = 0;
  
  do {
    graph = (*Graphmat)(newEstimate.theta);
    oldEstimate = newEstimate;
    newEstimate = mStep(y, oldEstimate, graph, wMatrix(y, oldEstimate), lambda);
    
  } while (counter++ < maxRep && oldEstimate.distance(newEstimate) >= epsilon);
  
  if (verbose) {
    Rcpp::Rcout << "Total MEM iterations: " << counter << ".\n";
  }
  
  return newEstimate;
}

// Creates a matrix of w_ik values.
MatrixXd wMatrix(const MatrixXd& y, const Psi& psi) {
  MatrixXd result(n, K);
  double acc;

  for (int i = 0; i < n; i++) {
    acc = 0.0;

    for (int k = 0; k < K; k++) {
      result(i, k) = psi.pii(k) * (*density)(y.row(i), psi.theta.col(k), psi.sigma);
      acc += result(i, k);
    }

    result.row(i) /= acc;
  }

  return result;
}

// M-Step of the Modified EM Algorithm.
Psi mStep(const MatrixXd& y, const Psi& psi, const MatrixXd& graph,  const MatrixXd& wMtx, double lambda) {
  Psi result;
  int i, k;
  result.pii = VectorXd::Zero(K);
  
  if(!arbSigma) {
    double acc;
    
    for (k = 0; k < K; k++) {
      acc = 0.0;
      
      for (i = 0; i < n; i++) {
        acc += wMtx(i, k);
      }
      
      result.pii(k) = (acc + ck) / (n + K * ck);
    }
    
    result.sigma = psi.sigma;
    
    // The following is currently hard-coded for updating 
    // the common structure parameter ("sigma")
    // in multivariate location Normal mixtures.
    // The notation is that of McLachlan and Peel 
    // (Finite Mixture Models, Chapter 3). 
  } else {   
    double Tk1Sum;      
    MatrixXd Tk2Sum(N, 1);
    MatrixXd Tk3Sum(N, N);
    
    result.sigma = MatrixXd::Zero(N, N);
    
    for (k = 0; k < K; k++) {
      Tk1Sum = wMtx.col(k).sum();
      Tk2Sum = MatrixXd::Zero(N, 1); //(wMtx.col(k).transpose() * y).transpose();
      Tk3Sum = MatrixXd::Zero(N, N);
      
      result.pii(k) = (Tk1Sum + ck) / (n + K * ck);
      
      for (i = 0; i < n; i++) {
        Tk2Sum += wMtx(i, k) * y.row(i).transpose();
        Tk3Sum += wMtx(i, k) * /*yOuterProd[i]; */ y.row(i).transpose() * y.row(i);
      }
      
      result.sigma += Tk3Sum - (1.0 / Tk1Sum) * Tk2Sum * Tk2Sum.transpose();
    }
    
    result.sigma /= n;
  }
  
  result.theta = psi.theta;
  
  result.theta = admm(y, result, graph, wMtx, lambda);
  
  return result;
}

//Revised ADMM algorithm
MatrixXd admm(const MatrixXd& y, const Psi& psi, const MatrixXd& graph, const MatrixXd& wMtx, double lambda){
  int k,j, counter = 0;
  double phi;
  MatrixXd newTheta = MatrixXd::Zero(D, K), 
    oldTheta(D, K), 
    Eta(D,K*K), 
    oldU = MatrixXd::Zero(D,K*K),
    newU(D,K*K);
  
  //Initialize theta
  
  oldTheta = psi.theta;
  for(k = 0; k < K; k++) {
    for(j = 0; j < K; j++){
      Eta.col(k+K*j) = psi.theta.col(k);
    }}
  
  do {
    oldTheta = newTheta ; 
    newTheta = (*updateTheta)(y, psi.sigma, graph, wMtx, Eta, oldU);
    
    
    //update Eta
    for(k = 0; k < K; k++) {
      for(j = 0; j < K; j++){
        if (graph(k,j)==1){
          phi = etamax(newTheta.col(k)+oldU.col(k+K*j)-newTheta.col(j)-oldU.col(j+K*k), lambda);
          Eta.col(k+K*j) = phi*(newTheta.col(k)+oldU.col(k+K*j))+(1-phi)*(newTheta.col(j)+oldU.col(j+K*k));
        }}}
    
    //update U
    for(k = 0; k < K; k++) {
      for(j = 0; j < K; j++){
        if (graph(k,j)==1){
          newU.col(k+K*j)=oldU.col(k+K*j)+newTheta.col(k)-Eta.col(k+K*j);
        }}}
    oldU=newU;
    
  } while (counter++ < maxadmm ||(oldTheta - newTheta).norm() < delta);
  return newTheta;
}

// Proximal Gradient Descent Algorithm. 
MatrixXd pgd(const MatrixXd& y, const Psi& psi, const MatrixXd& wMtx, double lambda){
  int k, counter = 0, restartCounter;
  double u;
  MatrixXd newEta = MatrixXd::Zero(D, K), 
           newTheta(D, K), 
           oldEta(D, K), 
           wMatSum, 
           z(D, K), 
           incrMatrix(D, K),
           newEtaTransf;
  
  VectorXd wMtxSums(K);
  
  // Initialize eta
  oldEta = psi.theta * invEtaTransf;

  u = uBound;

  newTheta = psi.theta;      

  wMatSum    = wMtx.transpose() * ones;
  incrMatrix = transformedData * wMtx * upTransf.transpose() - 
               (*gradB)(newTheta, psi.sigma) * wMatSum.triangularView<Lower>();

  for(k = 0; k < K; k++) {
    wMtxSums(k) = wMtx.col(k).sum();
  }

  do {
    Rcpp::checkUserInterrupt();
    restartCounter = 0.0;

    z = oldEta + u * incrMatrix;

    newEta.col(0) = z.col(0);

    // SCAD or MCP. 
    if (lambda > 0 && (penalty == 1 || penalty == 2)) {
      for (k = 1; k < K; k++) {
        newEta.col(k) = (*updateEta)(u, z.col(k), lambda, a);
      }

    // LLA SCAD or MCP.
    } else if (lambda > 0 && (penalty == 4 || penalty == 5)) {
      for (k = 1; k < K; k++) {
        newEta.col(k) = (*updateEtaLLA)(u, z.col(k), oldEta.col(k), lambda, a);
      }

    // Adaptive Lasso.
    } else if (lambda > 0 && penalty == 3) {
      for (k = 1; k < K; k++) {
        newEta.col(k) = (*updateEta)(u, z.col(k), lambda, adaptiveLassoWeights(k));
      }

    } else {
      newEta = z;
    }

    if(!(*constrCheck)(newEta * upTransf) || !uCheck(u, y, oldEta, newEta, psi.sigma, wMtxSums)) {
      if(restartCounter++ > maxPgd) break;

      u *= 0.5;
      continue;

    } else {
      if (2 * u < uBound) {
        u *= 2;

      } else {
        u = uBound;
      }

      newTheta = newEta * upTransf;

      if((oldEta - newEta).norm() < delta) {
         break;
      }

      oldEta     = newEta;
      incrMatrix = transformedData * wMtx * upTransf.transpose() - 
                   (*gradB)(newTheta, psi.sigma) * wMatSum.triangularView<Lower>();
    }

  } while(counter++ < maxPgd);

  return newTheta;
}

// Learning the tuning parameter u.
bool uCheck(double u, const MatrixXd& y, const MatrixXd& oldEta, const MatrixXd& newEta, const MatrixXd& sigma, const VectorXd& wMtxSums) {
  double t1 = 0.0, t2 = 0.0;
  
  MatrixXd grad(D, K);

  MatrixXd oldTheta = oldEta * upTransf,
           newTheta = newEta * upTransf; 

  grad = (*gradB)(oldTheta, sigma);

  for (int k = 0; k < K; k++){ 
    t2 += wMtxSums(k) * (b(newTheta.col(k), sigma) - b(oldTheta.col(k), sigma) - (grad.col(k).transpose() * (newTheta.col(k) - oldTheta.col(k)))(0, 0)); 
  }

  t1 = (1 / (2 * u)) * (newEta - oldEta).squaredNorm();

  return (t1 >= t2);
}

// Log-likelihood function. 
double fullLogLikFunction(const MatrixXd& y, const MatrixXd& theta, const VectorXd& pii, const MatrixXd& sigma){
  double temp, loglikSum = 0.0;

  for (int i = 0; i < n; i++) {
    temp = 0.0;

    for (int k = 0; k < K; k++) {
      temp += pii(k) * (*density)(y.row(i), theta.col(k), sigma);
    }

    loglikSum += log(temp);
  }

  return loglikSum;
}

double logLikFunction(const MatrixXd& y, const Psi& psi){
  return fullLogLikFunction(y, psi.theta, psi.pii, psi.sigma);
}

Rcpp::List estimateSequence(const MatrixXd& y, const Psi& startingVals, const VectorXd& lambdaList){
  Psi psi = startingVals, minPsi;
  int i, k;
  MatrixXd transfTheta;
  
  Rcpp::List estimates;
  Rcpp::NumericVector rbicVals, orders, loglikVals;

  if(penalty == 3) {
    adaptiveLassoWeights = VectorXd::Zero(K);
    adaptiveLassoWeights(0) = 1; //1/(psi.theta.col(0).norm() + 0.0001);
    psi = mem(y, psi, alStart);
    for (k = 1; k < K; k++) {
      adaptiveLassoWeights(k) = 1/(0.0001 + (psi.theta.col(k) - 
                                psi.theta.col(k-1)).norm());
    }
  }

  for (i = 0; i < lambdaList.size(); i++) {
    Rcpp::List thisEstimate;
    Rcpp::List th(K);
    Rcpp::NumericVector pii;
    Rcpp::CharacterVector names(K);

    if (verbose) 
      Rcpp::Rcout << "Lambda " << lambdaList(i) << ".\n";

    try {
 
      psi = mem(y, psi, lambdaScale * lambdaList(i));

      if (verbose) 
        Rcpp::Rcout << "Estimate: \n" << invTransf(psi.theta, psi.sigma) << "\n\n";
   
    } catch (const char* error) {
      throw error;
    } 
   
    pii = Rcpp::wrap(psi.pii);

    transfTheta = invTransf(psi.theta, psi.sigma);

    for(k = 0; k < K; k++) {
      std::ostringstream oss1;
      oss1 << k + 1;
      names[k] = "th" + oss1.str();
      th[k]    = transfTheta.col(k);
    }

    Rcpp::DataFrame theta(th);
    theta.attr("names") = names;

    std::ostringstream oss;
    oss << lambdaList(i);

    if (i == 0) {
      thisEstimate["ck"] = ck;
    }

    thisEstimate["lambda"] = lambdaList(i); 
    thisEstimate["order"]  = frequency(psi.theta);
    thisEstimate["pii"]    = pii;

    switch (modelIndex) {
      case 1: thisEstimate["mu"]    = transfTheta;
              thisEstimate["sigma"] = psi.sigma;
              break;
        
      case 2: thisEstimate["mu"]    = transfTheta.row(0);
              thisEstimate["sigma"] = transfTheta.row(1);
              break;

      default: thisEstimate["theta"] = transfTheta;
    }
    
    estimates.push_back(thisEstimate);
  }

  return estimates;
}
