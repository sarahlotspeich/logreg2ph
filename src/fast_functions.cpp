#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
MatrixXd fastMatrixMultiply(const MatrixXd& mat1, const MatrixXd& mat2)
{
	return mat1 * mat2;
}

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


// [[Rcpp::export]]
MatrixXd elementMultiply(const MatrixXd &mat, const VectorXd& vec)
{
	// int matSize = mat.rows() * mat.cols();
	// int vecSize = vec.size();
	// if (matSize != vecSize)
	// {
	// 	if (matSize % vecSize == 0)
	// 	{
	// 		VectorXd longVec(matSize);
	// 		for (int i = 0; i < matSize / vecSize; ++i)
	// 		{
	// 			longVec << vec;
	// 		}
	// 		vec = longVec;
	// 	}
	// }

	const MatrixXd prod = mat.cwiseProduct(vec);
	NumericMatrix m = eigenToRcpp(prod);
	m.attr("dim") = Dimension(mat.rows(), mat.cols());
	return as<MatrixXd>(m);
}

// [[Rcpp::export]]
VectorXd elementMultiplyVec(const VectorXd &mat, const VectorXd &vec)
{
	// Identical to mat.array() * vec.array()
	return mat.cwiseProduct(vec);
}



// [[Rcpp::export]]
VectorXd calculateMu(
	const MatrixXd& design_mat,
	const MatrixXd& prev)
{
	MatrixXd mu = design_mat * prev * -1;
	VectorXd mu1 = mu.array().exp();
	NumericVector muExp = eigenToRcpp(mu1);
	const VectorXd muVector = as<VectorXd>(NumericVector(muExp / (1 + muExp)));

	return muVector;
}

// [[Rcpp::export]]
NumericVector calculateGradient(
	const VectorXd& w_t_original,
	const int& n,
	const NumericMatrix& design_mat_R,
	const NumericVector& Y_col_R,
	const VectorXd& muVector,
	const bool& modifyW_T = true)
{
	// Convert from Rcpp to RcppEigen for some faster calculations
	const MatrixXd design_mat = as<MatrixXd>(design_mat_R);
	const VectorXd Y_col = as<VectorXd>(Y_col_R);

	// Put n 1's in front of w_t
	VectorXd w_t(n + w_t_original.size());
	if (modifyW_T)
	{
		w_t << VectorXd::Constant(n, 1), w_t_original;
	}
	else
	{
		VectorXd w_t = w_t_original;
	}

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
	const VectorXd& w_t_original,
	const NumericVector& muVector,
	const int & n,
	const bool& modifyW_T = true)
{
	VectorXd w_t(n + w_t_original.size());
	if (modifyW_T)
	{
		w_t << VectorXd::Constant(n, 1), w_t_original;
	}
	else
	{
		VectorXd w_t = w_t_original;
	}

	// post_multiply <- c(w_t * muVector * (muVector - 1)) * gamma_design_mat
	const NumericVector mus = muVector * (muVector - 1);
	const NumericVector c = mus * eigenToRcpp(w_t);
	const NumericMatrix post_multiply = matrixTimesVector(design_mat, c);

	return applyHessian(post_multiply, design_mat);
}

