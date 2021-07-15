#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <stdlib.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace Eigen;
using namespace std;

typedef Map<MatrixXd> MapMatrix;
typedef Map<VectorXd> MapVector;


// CONVERTING TYPE HELPER FUNCTIONS
// Since these work at the memory level, these are lossless

NumericMatrix eigenToRcpp(const MatrixXd& mat)
{
  SEXP s = wrap(mat);
  NumericMatrix ret(s);
  return ret;
}

NumericVector eigenToRcpp(const VectorXd& vec)
{
  SEXP s = wrap(vec);
  NumericVector ret(s);
  return ret;
}

Eigen::MatrixXd armaToEigen(arma::mat arma_A) {

  Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(),
    arma_A.n_rows,
    arma_A.n_cols);

  return eigen_B;
}

arma::mat eigenToArma(MatrixXd eigen_A)
{
  arma::mat arma_B = arma::mat(eigen_A.data(),
    eigen_A.rows(),
    eigen_A.cols(),
    true,
    false);
  return arma_B;
}

arma::mat rcppToArma(NumericMatrix mat)
{
  arma::mat armaMat = arma::mat(mat.begin(),
    mat.nrow(),
    mat.ncol(),
    false);
  return armaMat;
}

NumericMatrix armaToRcpp(arma::mat arma_A)
{
  return as<NumericMatrix>(wrap(arma_A));
}

NumericVector armaToRcpp(arma::vec arma_A)
{
  return as<NumericVector>(wrap(arma_A));
}


// RECREATING R FUNCTIONALITY

// [[Rcpp::export]]
MatrixXd fastMatrixMultiply(const MapMatrix& mat1, const MapMatrix& mat2)
{
  return mat1 * mat2;
}

// Don't export this overloaded function, for use with non-mapped matrices only
MatrixXd fastMatrixMultiply(const MatrixXd& mat1, const MatrixXd& mat2)
{
  return mat1 * mat2;
}

// R's built-in behavior of multiplying a matrix by a vector (element-wise)
// [[Rcpp::export]]
NumericMatrix matrixTimesVector(NumericMatrix mat, NumericVector vec){
  NumericMatrix res(mat.rows(), mat.cols());
  int index = 0;
  for (int i = 0; i < mat.cols(); ++i){
    for (int j = 0; j < mat.rows(); ++j){
      res( j , i ) = mat( j , i ) * vec[index++ % vec.size()];
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::mat matTimesVec(arma::mat mat, arma::vec v) {
// Ensure the vector is the right length
  if (v.n_elem < mat.n_rows)
  {
    arma::vec oldV = v;
    for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
    {
      v = join_vert(v, oldV);
    }
  }

// Multiply each col by the vector
  mat.each_col() %= v;
  if (mat.has_nan())
  {
    Rcout << "matTimesVec IS THE PROBLEM";
  }
  return mat;
}

// R's rep() function
NumericVector repeat(NumericVector x, NumericVector each) {

  NumericVector myvector(sum(each));

  for (int i = 0; i < each.size(); ++i)
  {
    int ind = 0;
    for (int j = 0; j < each[i]; ++j)
    {
      myvector[ind++] = x[i];
    }
  }
  return myvector;
}

// [[Rcpp::export]]
MatrixXd elementMultiply(const MapMatrix &mat, const MapVector& vec)
{
  const MatrixXd prod = mat.cwiseProduct(vec);
  NumericMatrix m = eigenToRcpp(prod);
  m.attr("dim") = Dimension(mat.rows(), mat.cols());
  return as<MatrixXd>(m);
}

// [[Rcpp::export]]
VectorXd elementMultiplyVec(const MapVector& mat, const MapVector& vec)
{
// Identical to mat.array() * vec.array()
  return mat.cwiseProduct(vec);
}


// TRANSLATING PACKAGE FUNCTIONS TO CPP FOR SPEED BOOST

// [[Rcpp::export]]
arma::vec lengthenWT(
  const arma::vec& w_t_original,
  const int& n,
  const bool& modifyW_T = true)
{
  if (!modifyW_T)
    return w_t_original;

// Put n 1's in front of w_t
  arma::vec w_t(n + w_t_original.n_elem);

  w_t = join_vert(arma::ones<arma::vec>(n), w_t_original);


  return w_t;
}

// [[Rcpp::export]]
arma::vec calculateMu(
  const arma::mat& design_mat,
  const arma::mat& prev)
{
  arma::mat mu = (design_mat * prev) * -1;
  arma::vec mu1 = exp(mu).as_col();
  return mu1 / (1 + mu1);
}

// [[Rcpp::export]]
arma::vec calculateGradient(
  const arma::vec& w_t_original,
  const int& n,
  const arma::mat& design_mat,
  const arma::vec& Y_col,
  const arma::vec& muVector,
  const bool& modifyW_T = false)
{
// Put n 1's in front of w_t
  arma::vec w_t = lengthenWT(w_t_original, n, modifyW_T);

// Calculate gradient
  arma::vec sumsVector = Y_col - 1 + muVector;

// Convert to NumericMatrix for Rcpp sugar's colSums
  arma::mat temp = matTimesVec(design_mat, sumsVector);
  temp = matTimesVec(temp, w_t);

// Sum's default behavior on matrices is like colSums
  arma::rowvec gradient = sum(temp);


  return reshape(gradient, gradient.n_elem, 1);

}

// [[Rcpp::export]]
arma::mat applyHessian(
  const arma::mat& post_multiply,
  const arma::mat& design_mat_R)
{
// Make a square matrix with size design_mat_R.cols()
  arma::mat hessian_gamma(design_mat_R.n_cols, design_mat_R.n_cols);

// For each col of design_mat_R, we multiply that column by a col from post_multiply then take the sum
// We end up with a matrix of every combo of columns
  for (int i = 0; i < design_mat_R.n_cols; ++i)
  {
    for (int j = 0; j < design_mat_R.n_cols; ++j)
    {
// Get references to the columns
      const arma::vec designCol = design_mat_R.col(i);
      const arma::vec thisCol = post_multiply.col(j);

// Multiply columns
      const arma::vec productVector = thisCol % designCol;

// Include this sum in our matrix
      hessian_gamma(i, j) = sum(productVector);
    }
  }

  return hessian_gamma;
}

// [[Rcpp::export]]
arma::mat calculateHessian(
  const arma::mat& design_mat,
  const arma::vec& w_t_original,
  const arma::vec& muVector,
  const int & n,
  const bool& modifyW_T = false)
{
  arma::vec w_t = lengthenWT(w_t_original, n, modifyW_T);

// post_multiply = c(w_t * muVector * (muVector - 1)) * gamma_design_mat
  const arma::vec mus = muVector % (muVector - 1);
  const arma::vec c = mus % w_t;
  const arma::mat post_multiply = matTimesVec(design_mat, c);

  return applyHessian(post_multiply, design_mat);
}

// VectorXd sumOverLogPTheta(
// 	const MatrixXd& comp_dat_theta_pred,
// 	const MatrixXd& comp_dat_Y_val,
// 	const VectorXd& theta)
// {
// 	// cbind(int=1, comp_dat_theta_pred)
// 	// Append a col of 1's to the left of our matrix
// 	DataFrame pY_X = DataFrame::create(_["int"] = eigenToRcpp(VectorXd::Constant(1, comp_dat_theta_pred.rows())),
// 		eigenToRcpp(comp_dat_theta_pred * theta));
// 	pY_X *= -1;
// 	pY_X = (pY_X.array().exp() + 1) ^ -1;

// 	return pY_X;
// }

// [[Rcpp::export]]
arma::vec pYstarCalc(
  const bool& errorsY,
  const arma::mat& gamma_design_mat,
  const int& n,
  const arma::mat& prev_gamma,
  const arma::mat& comp_dat_all,
  const int& Y_unval_index  )
{
  arma::vec pYstar;
  arma::mat mu_gamma;
  if (errorsY)
  {
// same as gamma_design_mat[-c(1:n),]
// get the elements of gamma_design_mat excluding the first n rows
    arma::mat filtered_gamma_design_mat = gamma_design_mat.rows(n, gamma_design_mat.n_rows-1);

    mu_gamma = filtered_gamma_design_mat * prev_gamma;
    pYstar = 1 / (1 + exp(mu_gamma * -1));

    arma::vec checkVector = comp_dat_all.col(Y_unval_index).rows(n, comp_dat_all.n_rows-1);
    for (int i = 0; i < pYstar.size(); ++i)
    {
      if (checkVector(i) == 0)
        pYstar(i) = 1 - pYstar(i);
    }

  }

  return pYstar;
}

// [[Rcpp::export]]
arma::mat pXCalc(
  const int& n,
  const arma::mat& comp_dat_all,
  const bool& errorsX,
  const bool& errorsY,
  const arma::mat& prev_p,
  const arma::uvec& indices,
  const arma::uvec& Bspline_index)
{
  arma::mat pX(1,1);
  IntegerVector rowIndices = seq(n, comp_dat_all.n_rows-1);
  if (errorsX and errorsY) {
// p_kj
// need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
// multiply by the B-spline terms

// these indices need to be repeated 2x
    arma::mat joinedPrevP = join_vert(prev_p.rows(indices), prev_p.rows(indices));
    arma::mat compDatSubmat = comp_dat_all.submat(as<arma::uvec>(rowIndices), Bspline_index);

// Rcout << "joinedPrevP: " << joinedPrevP.n_rows << " x " << joinedPrevP.n_cols << " | compDatSubmat: " << compDatSubmat.n_rows << " x " << compDatSubmat.n_cols;

// element-wise multiplication
    pX =  joinedPrevP % compDatSubmat;
//  p_kj
  } else if (errorsX) {
// p_kj
// need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
// multiply by the B-spline terms
    pX = prev_p.rows(indices) % comp_dat_all.submat(as<arma::uvec>(rowIndices), Bspline_index);
//  p_kj
  }

  return pX;
}

// [[Rcpp::export]]
List conditionalExpectations(const bool& errorsX,
  const bool& errorsY,
  const arma::mat& pX,
  const arma::vec& pY_X,
  const arma::vec& pYstar,
  const arma::uvec& psiDenomIndex,
  const int& nDiff,
  const int& m)
{
  arma::vec w_t(1);
  arma::mat u_t(1,1);
  arma::mat psi_num;
  arma::mat psi_t;
  if (errorsY and errorsX)
  {

// P(Y|X,C)P(Y*|X*,Y,X,C)p_kjB(X*)
    arma::mat firstPart = matTimesVec(pX, pY_X);
    psi_num = matTimesVec(firstPart, pYstar);

// Update denominator
// Sum up all rows per id (e.g. sum over xk/y)
/*
psi_denom looks like: (for N-n = 750)
psi_num[0,] + psi_num[750,] + psi_num[1500,] + ...
psi_num[1,] + psi_num[751,] + ...
...
psi_num[749,] + psi_num[1499] + ...
*/
    arma::mat combined = join_vert(psi_num.rows(psiDenomIndex), psi_num.rows(psiDenomIndex)) ;
    arma::mat psi_denom(nDiff, combined.n_cols, arma::fill::zeros);

    for (int i = 0; i < psi_num.n_rows; ++i)
    {
      psi_denom.row(i % nDiff) += psi_num.row(i);
    }

// psi_denom is correct here July 8

// Then sum over the sn splines
    psi_denom = sum(psi_denom, 1);

// Avoid NaN resulting from dividing by 0
    for (int i = 0; i < psi_denom.n_rows; ++i)
    {
      if (psi_denom[i] == 0)
        psi_denom[i] = 1;
    }



// And divide them!
// We need to raise psi_denom to the power of -1 to essentially divide
// Otherwise, the matrix sizes are incompatible for built-in division
    arma::mat divisor = pow(psi_denom, -1);
    psi_t = matTimesVec(psi_num, divisor);

// Update the w_kyi for unvalidated subjects
// by summing across the splines/ columns of psi_t
// Equivalent to rowSums()
    w_t = sum(psi_t, 1);

// Rcout << "w_t is correct" << endl;

// Update the u_kji for unvalidated subjects
// by summing over Y = 0/1 w/i each i, k
// add top half of psi_t (y = 0) to bottom half (y = 1)

    u_t.reshape(size(psi_t));

    try
    {
      u_t = psi_t.head_rows(m * nDiff) + psi_t.tail_rows(m * nDiff);
    }
    catch(std::exception e)
    {
      Rcout << "Something has gone wrong!" << endl;
      return List::create(Named("u_t") = u_t, Named("w_t") = w_t, Named("psi_t")=psi_t);
    }

// Rcout << "all done" <<endl;
  }
  else if (errorsX)
  {
    Rcout << "JUST ERRORS X!";
    // P(Y|X,C)p_kjB(X*)
    psi_num = matTimesVec(pX, pY_X);
    // Update denominator
    // Sum up all rows per id (e.g. sum over xk)
    arma::mat psi_denom(nDiff, psi_num.n_cols, arma::fill::zeros);
    for (int i = 0; i < psi_num.n_rows; ++i)
    {
      psi_denom.row(i % nDiff) += psi_num.row(i);
    }
    // Then sum over the sn splines
    psi_denom = sum(psi_denom, 1);
    // Avoid NaN resulting from dividing by 0
    for (int i = 0; i < psi_denom.n_rows; ++i)
    {
      if (psi_denom[i] == 0)
        psi_denom[i] = 1;
    }
    // And divide them!
    psi_t = matTimesVec(psi_num, pow(psi_denom, -1));

    // Update the w_kyi for unvalidated subjects
    // by summing across the splines/ columns of psi_t
    w_t = sum(psi_t, 1);
  }
  else if (errorsY)
  {
    Rcout << "JUST ERRORS Y!";
    // P(Y|X,C)P(Y*|Y,X,C)
    arma::colvec psi_num = pY_X * pYstar;
    // Sum up all rows per id (e.g. sum over y)
    arma::mat psi_denom(size(psi_num), arma::fill::zeros);

    for (int i = 0; i < psi_num.n_rows; ++i)
    {
      psi_denom.row(i % nDiff) += psi_num.row(i);
    }
    // Avoid NaN resulting from dividing by 0
    for (int i = 0; i < psi_denom.n_rows; ++i)
    {
      if (psi_denom[i] == 0)
        psi_denom[i] = 1;
    }
    // And divide them!
    psi_t = matTimesVec(psi_num, pow(psi_denom, -1));
    // Update the w_kyi for unvalidated subjects
    w_t = psi_t;
  }

  return List::create(Named("w_t")=w_t, Named("u_t")=u_t, Named("psi_t")=psi_t);
}

// [[Rcpp::export]]
List profileOutLoop(
  const int& MAX_ITER,
  const bool& errorsX,
  const bool& errorsY,
  const int& Y_unval_index,
  const arma::uvec& Bspline_index,
  const int& N,
  const int& n,
  const int& m,
  const arma::mat& gamma_design_mat,
  const arma::mat& prev_gamma_R,
  const arma::mat& comp_dat_all,
  const arma::mat& prev_p_R,
  const arma::uvec& indices,
  const arma::vec& pY_X,
  const arma::uvec& psiDenomIndex,
  const arma::mat& gamma0,
  const arma::mat& p0,
  const float& TOL,
  const arma::mat& p_val_num)
{
  int it = 1;
  bool CONVERGED = false;

  arma::mat new_gamma;
  arma::mat new_p;
  arma::mat psi_t;

// If we don't create deep copies of these, then the cpp modifies them, messing up R
  arma::mat prev_gamma(prev_gamma_R);
  arma::mat prev_p(prev_p_R);

// Estimate gamma/p using EM
  while(it <= MAX_ITER && !CONVERGED)
  {
// E Step
// P(Y*|X*,Y,X)
// Rcout << "start loop" << endl;

    arma::vec pYstar = pYstarCalc(errorsY, gamma_design_mat, n, prev_gamma, comp_dat_all, Y_unval_index);
// pYstar is correct at this point - June 30

// Rcout << "pYstar is correct" << endl;
// P(Y*|X*,Y,X)
// ###################################################################
// P(X|X*)

    arma::mat pX = pXCalc(n, comp_dat_all, errorsX, errorsY, prev_p, indices, Bspline_index);
// pX is exactly correct at this point - July 6
// Rcout << "pX is correct" << endl;


//   P(X|X*)
//   Estimate conditional expectations
    List temp = conditionalExpectations(errorsX, errorsY, pX, pY_X, pYstar, psiDenomIndex, N-n, m);
    arma::vec w_t = temp["w_t"];
    arma::mat u_t = temp["u_t"];
    arma::mat psi_t = temp["psi_t"];

// u_t, w_t are correct here July 12

//  Estimate conditional expectations
//  E Step
//   ###################################################################

//   ###################################################################
// M Step
//   ###################################################################
    bool gamma_conv = false;
    bool p_conv = false;

    if (errorsY)
    {
// Update gamma using weighted logistic regression
      w_t = lengthenWT(w_t, n);
      const arma::vec muVector = calculateMu(gamma_design_mat, prev_gamma);
      const arma::mat gradient_gamma = calculateGradient(w_t, n, gamma_design_mat, comp_dat_all.col(Y_unval_index), muVector);
      arma::mat hessian_gamma = calculateHessian(gamma_design_mat, w_t, muVector, n);

// Rcout << "hessian_gamma" << endl;

// return List::create(w_t, muVector, gradient_gamma, hessian_gamma);
// Gradient
// Hessian
      try
      {
// If missing, b is taken to be an identity matrix and solve will return the inverse of a.
// return List::create(Named("prev_gamma")=prev_gamma,Named("hessian_gamma")=(hessian_gamma), Named("hessian_gamma.i()")=inv(hessian_gamma), Named("gradient_gamma")=gradient_gamma);
        new_gamma = prev_gamma - (inv(hessian_gamma) * gradient_gamma);
// Rcout << it << ": " << new_gamma << endl;
      }
      catch(std::exception e)
      {
        Rcout << "exception thrown" << endl;
        new_gamma.reshape(0,0);
      }

      if (new_gamma.has_nan())
      {
        Rcout << "we need glm!!" << endl;
        return List::create(new_gamma);
// new_gamma = matrix( glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1);
      }
// Check for convergence
      gamma_conv = abs(new_gamma - prev_gamma).max() < TOL;

//  Update gamma using weighted logistic regression
    } else { gamma_conv = TRUE; }

// Update {p_kj}
    if (errorsX)
    {
      arma::mat sumMatrix(m, u_t.n_cols, arma::fill::zeros);
      if (errorsY)
      {
// Update numerators by summing u_t over i = 1, ..., N
        for (int i = 0; i < m; ++i)
          for (int j = 0; j < N-n; ++j)
            sumMatrix.row(i) += u_t.row(i*m + N-n);

        }
        else
        {

// Update numerators by summing u_t over i = 1, ..., N
// Difference between here and above is u_t vs psi_t
          sumMatrix.set_size(m, psi_t.n_cols);
          sumMatrix.fill(0);
          for (int i = 0; i < m; ++i)
            for (int j = 0; j < N-n; ++j)
              sumMatrix.row(i) += psi_t.row(i*m + N-n);

          }

          sumMatrix = sort(sumMatrix);
          arma::mat new_p_num = p_val_num + sumMatrix;
          new_p = matTimesVec( new_p_num.t(), pow(sum(new_p_num).t(), -1) ).t();

// Check for convergence
          p_conv = abs(new_p - prev_p).max() < TOL;

        }
        else { p_conv = TRUE; }
// Rcout << "new_p" << endl;

//  Update {p_kj}
// ###################################################################
// M Step
// ###################################################################

        CONVERGED = gamma_conv && p_conv;

// Update values for next iteration
        ++it;
        if (errorsY) { prev_gamma = new_gamma; }
        if (errorsX) { prev_p = new_p; }
//   Update values for next iteration
      }
// Rcout << "loop all done" << endl;
// Rcout << size(new_gamma) << " " << size(new_p) << endl;
// return List::create(Named("psi_t")=psi_t,
//   Named("it")=it,
//   Named("CONVERGED")=CONVERGED);
      string CONVERGED_MSG = "";

      if(it == MAX_ITER && !CONVERGED)
      {
        CONVERGED_MSG = "MAX_ITER reached";
        if (errorsY)
        {
          new_gamma.reshape(gamma0.n_rows, 1);
          new_gamma.fill(NULL);
        }
        else { new_gamma.reshape(0,0); }

        if (errorsX)
        {
          new_p.reshape(p0.n_rows, p0.n_cols);
          new_p.fill(NULL);
        }
        else { new_p.reshape(0,0); }
      }

      if(CONVERGED)
        CONVERGED_MSG = "converged";
      if (!errorsY) { new_gamma.reshape(0,0); }
      if (!errorsX) { new_p.reshape(0,0); }

//  Estimate theta using EM
      return List::create(Named("psi_at_conv") = psi_t,
        Named("gamma_at_conv") = new_gamma,
        Named("p_at_conv") = new_p,
        Named("converged") = CONVERGED,
        Named("converged_msg") = CONVERGED_MSG);
    }


// [[Rcpp::export]]
    void testStuff(){

      arma::mat X ;
      arma::vec beta ;

      beta.resize ( 2 ) ;

      beta (0) = 1.0 ;
      beta (1) = 3.0 ;

      X.resize ( 3, 2 ) ;

      X (0,0) = 1.0 ;
      X (0,1) = 2.0 ;
      X (1,0) = 3.0 ;
      X (1,1) = 4.0 ;
      X (2,0) = 5.0 ;
      X (2,1) = 6.0 ;

      Rcout << X % beta << std::endl ;
    }
