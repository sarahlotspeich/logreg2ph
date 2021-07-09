#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace Eigen;
using namespace std;

typedef Map<MatrixXd> MapMatrix;
typedef Map<VectorXd> MapVector;


// CONVERTING TYPE HELPER FUNCTIONS

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
VectorXd lengthenWT(const MapVector& w_t_original,
	const int& n,
	const bool& modifyW_T = true)
{
	if (!modifyW_T)
		return w_t_original;

	// Put n 1's in front of w_t
	VectorXd w_t(n + w_t_original.size());
	w_t << VectorXd::Constant(n, 1), w_t_original;


	return w_t;
}

// [[Rcpp::export]]
VectorXd calculateMu(
	const MapMatrix& design_mat,
	const MapMatrix& prev)
{
	MatrixXd mu = design_mat * prev * -1;
	VectorXd mu1 = mu.array().exp();
	NumericVector muExp = eigenToRcpp(mu1);
	const VectorXd muVector = as<VectorXd>(NumericVector(muExp / (1 + muExp)));

	return muVector;
}

// [[Rcpp::export]]
NumericVector calculateGradient(
	const MapVector& w_t_original,
	const int& n,
	const NumericMatrix& design_mat_R,
	const NumericVector& Y_col_R,
	const MapVector& muVector,
	const bool& modifyW_T = false)
{
	// Convert from Rcpp to RcppEigen for some faster calculations
	const MatrixXd design_mat = as<MatrixXd>(design_mat_R);
	const VectorXd Y_col = as<VectorXd>(Y_col_R);

	// Put n 1's in front of w_t
	VectorXd w_t = lengthenWT(w_t_original, n, modifyW_T);

	// Define a constant vector to perform 'scalar' addition on vectors
	const VectorXd ones = VectorXd::Constant(Y_col.size(), 1);

	// Calculate gradient
	VectorXd sumsVector = Y_col - ones + muVector;
	NumericVector mulVector = eigenToRcpp(sumsVector);

	// Convert to NumericMatrix for Rcpp sugar's colSums
	NumericMatrix temp = matrixTimesVector(design_mat_R, mulVector);
	temp = matrixTimesVector(temp, eigenToRcpp(w_t));

	const NumericVector gradient = colSums(temp);
	return gradient;

}

// [[Rcpp::export]]
NumericMatrix applyHessian(
	const NumericMatrix& post_multiply,
	const NumericMatrix& design_mat_R)
{
	// Make a square matrix with size design_mat_R.cols()
	NumericMatrix hessian_gamma(design_mat_R.cols(), design_mat_R.cols());

  	// For each col of design_mat_R, we multiply that column by a col from post_multiply then take the sum
  	// We end up with a matrix of every combo of columns
	for (int i = 0; i < design_mat_R.cols(); ++i)
	{
		for (int j = 0; j < design_mat_R.cols(); ++j)
		{
      		// Get references to the columns
			const NumericVector designCol = design_mat_R(_, i);
			const NumericVector thisCol = post_multiply(_, j);

      		// Multiply columns
			const NumericVector productVector = thisCol * designCol;

      		// Sum this vector
			double sum = 0;
			for (int k = 0; k < productVector.size(); ++k)
				sum += productVector[k];

      		// Include this sum in our matrix
			hessian_gamma(i, j) = sum;
		}
	}

	return hessian_gamma;
}

// [[Rcpp::export]]
NumericMatrix calculateHessian(
	const NumericMatrix& design_mat,
	const MapVector& w_t_original,
	const NumericVector& muVector,
	const int & n,
	const bool& modifyW_T = false)
{
	VectorXd w_t = lengthenWT(w_t_original, n, modifyW_T);

	// post_multiply = c(w_t * muVector * (muVector - 1)) * gamma_design_mat
	const NumericVector mus = muVector * (muVector - 1);
	const NumericVector c = mus * eigenToRcpp(w_t);
	const NumericMatrix post_multiply = matrixTimesVector(design_mat, c);

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
arma::mat profileOutLoop(
  const int& it,
  const int& MAX_ITER,
  const bool& CONVERGED,
  const bool& errorsX,
  const bool& errorsY,
  const int& Y_unval_index,
  const arma::uvec& Bspline_index,
  const int& N,
  const int& n,
  const int& m,
  const arma::mat& gamma_design_mat,
  const MapMatrix& prev_gamma,
  const arma::mat& comp_dat_all,
  arma::mat& prev_p,
  const arma::uvec& indices,
  const NumericVector& pY_X,
  const arma::uvec& psiDenomIndex)
{
  // Estimate gamma/p using EM
	while(it <= MAX_ITER & !CONVERGED) {
  //   // E Step
  //   // P(Y*|X*,Y,X)
    Rcout << "start\n";
    arma::vec pYstar;
    Rcout << "initialize\n";
		if (errorsY)
    {
    	// same as gamma_design_mat[-c(1:n),]
    	// get the elements of gamma_design_mat excluding the first n rows
      MatrixXd filtered_gamma_design_mat = armaToEigen(gamma_design_mat.rows(n, gamma_design_mat.n_rows-1));

      MatrixXd mu_gamma = fastMatrixMultiply(filtered_gamma_design_mat, prev_gamma);
      pYstar = eigenToArma(1 / (1 + (mu_gamma * -1).array().exp()));

      VectorXd checkVector = armaToEigen(comp_dat_all.col(Y_unval_index).rows(n, comp_dat_all.n_rows-1));
      for (int i = 0; i < pYstar.size(); ++i)
      {
        if (checkVector(i) == 0)
          pYstar(i) = 1 - pYstar(i);
      }

    }

    // pYstar is correct at this point - June 30
Rcout << "pystar\n";

  //   //  P(Y*|X*,Y,X)
  //   ###################################################################
  //   // P(X|X*)
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
Rcout << "px" << endl;
        // pX is correct at this point - July 6

  //   P(X|X*)
  //   Estimate conditional expectations
    arma::vec w_t;
    arma::mat u_t;
    if (errorsY and errorsX) {
      // P(Y|X,C)P(Y*|X*,Y,X,C)p_kjB(X*)
      NumericMatrix firstPart = matrixTimesVector(armaToRcpp(pX),  pY_X);
      arma::mat psi_num(pX.n_rows, pX.n_cols);
      psi_num = rcppToArma(matrixTimesVector(firstPart, armaToRcpp(pYstar)));

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
      arma::mat psi_denom(N-n, combined.n_cols, arma::fill::zeros);

      for (int i = 0; i < N-n; ++i)
      {
        for (int j = i; j < psi_num.n_rows; j += N-n)
        {
          psi_denom.row(i) += psi_num.row(j);
        }
      }

      // psi_denom is correct here July 8

      // Then sum over the sn splines
      psi_denom = rowSums(armaToRcpp(psi_denom));

      // // Avoid NaN resulting from dividing by 0
      psi_denom(psi_denom == 0.0 - 1).fill(1);

      // And divide them!
      // We need to raise psi_denom to the power of -1 to essentially divide
      // Otherwise, the matrix sizes are incompatible
      arma::mat divisor = pow(psi_denom, -1);
      arma::mat psi_t = rcppToArma(matrixTimesVector(armaToRcpp(psi_num), armaToRcpp(divisor)));

      // // Update the w_kyi for unvalidated subjects
      // // by summing across the splines/ columns of psi_t
      w_t = rowSums(armaToRcpp(psi_t));
      // Update the u_kji for unvalidated subjects
      // by summing over Y = 0/1 w/i each i, k
      // add top half of psi_t (y = 0) to bottom half (y = 1)
      u_t.reshape(size(psi_t));
      u_t = psi_t.head_rows(m*(N-n)) + psi_t.tail_rows(m * (N - n));


    } else if (errorsX) {
      // P(Y|X,C)p_kjB(X*)
      // arma::mat psi_num = c(pY_X) % pX;
    //   // Update denominator
    //   // Sum up all rows per id (e.g. sum over xk)
      // arma::mat psi_denom = rowsum(psi_num, group = rep(seq(1, (N - n)), times = m));
    //   // Then sum over the sn splines
    //   psi_denom = rowSums(psi_denom);
    //   // Avoid NaN resulting from dividing by 0
    //   psi_denom[psi_denom == 0] = 1;
    //   // And divide them!
    //   arma::mat psi_t = psi_num / psi_denom;
    //   // Update the w_kyi for unvalidated subjects
    //   // by summing across the splines/ columns of psi_t
    //   w_t = rowSums(psi_t);
    // } else if (errorsY) {
    //   // P(Y|X,C)P(Y*|Y,X,C)
    //   psi_num = matrix(c(pY_X * pYstar), ncol = 1)
    //   // Sum up all rows per id (e.g. sum over y)
    //   psi_denom = rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2))
    //   // Avoid NaN resulting from dividing by 0
    //   psi_denom[psi_denom == 0] = 1
    //   // And divide them!
    //   psi_t = psi_num / psi_denom
    //   // Update the w_kyi for unvalidated subjects
    //   w_t = psi_t
    }

// u_t is correct here July 9
    return u_t;
  //   //  Estimate conditional expectations
  //   //  E Step
  //   ###################################################################

  //   ###################################################################
  //   // M Step
  //   ###################################################################
  //   if (errorsY) {
  //     // Update gamma using weighted logistic regression
  //     tic("M-step Y cpp")
  //     w_t = lengthenWT(w_t, n)
  //     muVector = calculateMu(gamma_design_mat, prev_gamma)
  //     gradient_gamma = calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[,c(Y_unval)], muVector)
  //     hessian_gamma = calculateHessian(gamma_design_mat, w_t, muVector, n)
  //     toc()


  //     // tic("M-step Y")
  //     // w_t = c(rep(1, n), w_t)
  //     // mu = gamma_design_mat %*% prev_gamma
  //     // gradient_gamma_R = matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(- mu) / (1 + exp(- mu)))) * gamma_design_mat)), ncol = 1)
  //     #
  //     //  Gradient
  //     // Hessian
  //     // post_multiply = c(w_t * muVector * (muVector - 1)) * gamma_design_mat
  //     // hessian_gamma_R = apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
  //     // if (max(abs(hessian_gamma - hessian_gamma_R)) > 1e-10)
  //     // {
  //     //   warning("cpp and R are significantly different!")
  //     //   browser()
  //     // }
  //     // toc()

  //     new_gamma = tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
  //                           error = function(err) {
  //                             matrix(NA, nrow = nrow(prev_gamma))
  //                           })
  //     if (any(is.na(new_gamma))) {
  //       suppressWarnings(new_gamma = matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
  //     }
  //     // Check for convergence
  //     gamma_conv = abs(new_gamma - prev_gamma) < TOL
  //     //  Update gamma using weighted logistic regression
  //   } else { gamma_conv = TRUE }
  //   ###################################################################
  //   // Update {p_kj}
  //   if (errorsX & errorsY) {
  //     tic("M-step both")
  //     // Update numerators by summing u_t over i = 1, ..., N
  //     new_p_num = p_val_num +
  //       rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
  //     new_p = t(t(new_p_num) / colSums(new_p_num))
  //     // Check for convergence
  //     p_conv = abs(new_p - prev_p) < TOL
  //     toc()
  //   } else if (errorsX) {
  //     tic("M-step X")
  //     // Update numerators by summing u_t over i = 1, ..., N
  //     new_p_num = p_val_num +
  //       rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
  //     new_p = t(t(new_p_num) / colSums(new_p_num))
  //     // Check for convergence
  //     p_conv = abs(new_p - prev_p) < TOL
  //     toc()
  //   }
  //   else { p_conv = TRUE }
  //   //  Update {p_kj}
  //   ###################################################################
  //   // M Step
  //   ###################################################################

  //   all_conv = c(gamma_conv, p_conv)
  //   if (mean(all_conv) == 1) { CONVERGED = TRUE }

  //   // Update values for next iteration
  //   it = it + 1
  //   if (errorsY) { prev_gamma = new_gamma }
  //   if (errorsX) { prev_p = new_p }
  //   //   Update values for next iteration
  }

  // tic("return profile_out")
  // if(it == MAX_ITER & !CONVERGED) {
  //   CONVERGED_MSG = "MAX_ITER reached"
  //   if (errorsY) { new_gamma = matrix(NA, nrow = nrow(gamma0), ncol = 1) } else { new_gamma = NA }
  //   if (errorsX) { new_p = matrix(NA, nrow = nrow(p0), ncol = ncol(p0)) } else { new_p = NA }
  // }
  // if(CONVERGED) CONVERGED_MSG = "converged"
  // if (!errorsY) { new_gamma = NA }
  // if (!errorsX) { new_p = NA }
  // toc()
  //  Estimate theta using EM
  // return(list("psi_at_conv" = psi_t,
  //             "gamma_at_conv" = new_gamma,
  //             "p_at_conv" = new_p,
  //             "converged" = CONVERGED,
  //             "converged_msg" = CONVERGED_MSG))
}
