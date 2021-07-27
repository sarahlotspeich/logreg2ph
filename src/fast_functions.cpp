// This definition allows us to do some big matrix multiplication (calculateHessian)
#define ARMA_64BIT_WORD 1


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
// Since these work at the memory level, these are lossless with time

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

// [[Rcpp::export]]
arma::mat matTimesVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  // Never an issue in this project, but if you're copy+pasting this code into yours, 
  // you may want to uncomment the error checking
  // if (v.n_elem < mat.n_rows)
  // {
  //   arma::vec oldV = v;
  //   for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
  //   {
  //     v = join_vert(v, oldV);
  //   }
  // }

  // Multiply each col by the vector
  mat.each_col() %= v;

  return mat;
}

// [[Rcpp::export]]
arma::mat matDivideVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  if (v.n_elem < mat.n_rows)
  {
    arma::vec oldV = v;
    for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
    {
      v = join_vert(v, oldV);
    }
  }

  // Divide each col by the vector
  mat.each_col() /= v;

  return mat;
}

// R's rep() function
NumericVector repeat(NumericVector x, NumericVector each)
{

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
  arma::mat retMat(design_mat.n_cols, design_mat.n_cols);
  retMat =  design_mat.t() * post_multiply;

  return retMat;
}

// [[Rcpp::export]]
arma::vec pYstarCalc(
  const bool& errorsY,
  const arma::mat& gamma_design_mat,
  const int& n,
  const int& excludeRows,
  const arma::mat& prev_gamma,
  const arma::mat& comp_dat_all,
  const int& Y_unval_index  )
{
  arma::vec pYstar;
  arma::mat mu_gamma;
  if (errorsY)
  {
    // same as gamma_design_mat[-c(1:n),]
    // get the elements of gamma_design_mat excluding the first excludeRows rows
    arma::mat filtered_gamma_design_mat = gamma_design_mat.rows(excludeRows, gamma_design_mat.n_rows-1);

    mu_gamma = filtered_gamma_design_mat * prev_gamma;
    pYstar = 1 / (1 + exp(mu_gamma * -1));

    arma::vec checkVector = comp_dat_all.col(Y_unval_index).rows(n, comp_dat_all.n_rows-1);
    for (int i = 0; i < pYstar.size(); ++i)
    {
      if (checkVector(i) == 0)
      {
        pYstar(i) = 1 - pYstar(i);
      }
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
  const arma::uvec& Bspline_index,
  const arma::uvec& rowIndices)
{
  arma::mat pX(1,1);
  arma::mat prevRows = prev_p.rows(indices);

  if (errorsX and errorsY) 
  {
    // need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    // multiply by the B-spline terms

    // these indices need to be repeated 2x
    arma::mat joinedPrevP = join_vert(prevRows, prevRows);

    // element-wise multiplication
    pX = joinedPrevP % comp_dat_all.submat(rowIndices, Bspline_index);
  } 
  else if (errorsX) 
  {
    // need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    // multiply by the B-spline terms
    pX = prevRows % comp_dat_all.submat(rowIndices, Bspline_index);
  }

  return pX;
}

// UNUSED
List conditionalExpectations(const bool& errorsX,
  const bool& errorsY,
  const arma::mat& pX,
  const arma::vec& pY_X,
  const arma::vec& pYstar,
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
    psi_num[749,] + psi_num[1499,] + ...
    */
    // arma::mat combined = join_vert(psi_num.rows(psiDenomIndex), psi_num.rows(psiDenomIndex)) ;
    arma::mat psi_denom(nDiff, psi_num.n_cols);

    for (int i = 0; i < psi_num.n_rows; ++i)
    {
      // This % is modulo division, not element-wise multiplication
      psi_denom.row(i % nDiff) += psi_num.row(i);
    }

    // psi_denom is correct here July 8

    // Then sum over the sn splines
    psi_denom = sum(psi_denom, 1);

    // Avoid NaN resulting from dividing by 0
    for (int i = 0; i < psi_denom.n_rows; ++i)
    {
      if (psi_denom[i] == 0)
      {
        psi_denom[i] = 1;
      }
    }



    // And divide them!
    psi_t = matDivideVec(psi_num, psi_denom);

    // Update the w_kyi for unvalidated subjects
    // by summing across the splines/ columns of psi_t
    // Equivalent to rowSums()
    w_t = sum(psi_t, 1);

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
      // Rcpp requires this exception handling, never encountered it
      Rcout << "Something has gone wrong while calculating u_t!" << endl;
      return List::create(Named("u_t") = u_t, Named("w_t") = w_t, Named("psi_t")=psi_t);
    }

    Rcout << "BOTH" << endl;

  }
  else if (errorsX)
  {


    // THIS BRANCH IS UNTESTED



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
    psi_t = matDivideVec(psi_num, psi_denom);

    // Update the w_kyi for unvalidated subjects
    // by summing across the splines/ columns of psi_t
    w_t = sum(psi_t, 1);
  }
  else if (errorsY)
  {


    // THIS BRANCH IS UNTESTED



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
    psi_t = matDivideVec(psi_num, psi_denom);

    // Update the w_kyi for unvalidated subjects
    w_t = psi_t;
  }

  return List::create(Named("w_t")=w_t, Named("u_t")=u_t, Named("psi_t")=psi_t);
}

// UNUSED
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
  const arma::mat& gamma0,
  const arma::mat& p0,
  const float& TOL,
  const arma::mat& p_val_num,
  const arma::uvec& prev_p_indices)
{
  int it = 1;
  bool CONVERGED = false;

  arma::mat new_gamma;
  arma::mat new_p;
  arma::mat psi_t(1,1);

// If we don't create deep copies of these, then the cpp modifies them, messing up R
  arma::mat prev_gamma(prev_gamma_R);
  arma::mat prev_p(prev_p_R);

// Estimate gamma/p using EM
  while(it <= MAX_ITER && !CONVERGED)
  {
    // E Step
    // P(Y*|X*,Y,X)

    arma::vec pYstar = pYstarCalc(errorsY, gamma_design_mat, n, n, prev_gamma, comp_dat_all, Y_unval_index);
    // pYstar is correct at this point - June 30


    // P(Y*|X*,Y,X)
    // ###################################################################
    // P(X|X*)

    arma::mat pX = pXCalc(n, comp_dat_all, errorsX, errorsY, prev_p, indices, Bspline_index, prev_p_indices);
    // pX is exactly correct at this point - July 6
    

    //   P(X|X*)
    //   Estimate conditional expectations
    List temp = conditionalExpectations(errorsX, errorsY, pX, pY_X, pYstar, N-n, m);
    arma::vec w_t = temp["w_t"];
    arma::mat u_t = temp["u_t"];
    arma::mat tt = temp["psi_t"];
    psi_t = tt;

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

    
      // Gradient
      // Hessian
      try
      {
        // If missing, b is taken to be an identity matrix and solve will return the inverse of a.
        new_gamma = prev_gamma - (inv(hessian_gamma) * gradient_gamma);
      }
      catch(std::exception e)
      {
        Rcout << "exception thrown" << endl;
        new_gamma.reshape(0,0);
      }

      if (new_gamma.has_nan())
      {
        // THIS BRANCH HAS NEVER BEEN ENTERED, SO NO IMPLEMENTATION OF GLM HAS BEEN MADE
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
        {
          sumMatrix.row(i) = sum(u_t.rows( i * (N - n), (i + 1) * (N - n) - 1 ));
        }
      }
      else
      {
        // Update numerators by summing u_t over i = 1, ..., N
        // Difference between here and above is u_t vs psi_t
        sumMatrix.set_size(m, psi_t.n_cols);
        sumMatrix.fill(0);
        for (int i = 0; i < m; ++i)
        {
          sumMatrix.row(i) = sum(psi_t.rows( i * (N - n), (i + 1) * (N - n) - 1 ));
        }

      }

      // sumMatrix = sort(sumMatrix);
      arma::mat new_p_num = p_val_num + sumMatrix;
      new_p = matDivideVec( new_p_num.t(), sum(new_p_num).t() ).t();

      // Check for convergence
      p_conv = abs(new_p - prev_p).max() < TOL;

      // new_p is correct here July 19
    }
    else { p_conv = TRUE; }

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

  if(CONVERGED) { CONVERGED_MSG = "converged"; }
  if (!errorsY) { new_gamma.reshape(0,0); }
  if (!errorsX) { new_p.reshape(0,0); }

  //  Estimate theta using EM
  return List::create(Named("psi_at_conv") = psi_t,
    Named("gamma_at_conv") = new_gamma,
    Named("p_at_conv") = new_p,
    Named("converged") = CONVERGED,
    Named("converged_msg") = CONVERGED_MSG);
}

