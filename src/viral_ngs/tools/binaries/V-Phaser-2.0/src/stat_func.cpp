//========================================================================
// Project     : VariantCaller
// Name        : stat_func.cpp
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================

#include "stat_func.h"

/* @brief	Given a prob vector [p], bonferroni factor, and  lamda value
 * 			of poisson distri, calculate threshold k such that f(X>=k)
 * 			is significant.
 * 			f(x) follows Poisson Binomial Distribution, which is calculated
 * 			either by an exact recursive function or poisson approximation
 */
int getThreshold (double lamda, const dvec_t& p, int bc) {

	int depth = p.size();
	int r = poisson_model(lamda, bc, depth);

	{
		//std::cout << "lamda = " << lamda << "\n";
		//std::cout << "vector (size: " << p.size() << ")\t";
		//debug_print_vec (p);
		//std::cout << "poisson threshold = " << r << "\n";
	}

	return r;

	/* recursive function: p[i] or q[i] cannot be 1 otherwise,
	 * std::log(1 + std::exp(newlog_p - newlog_q)) is buggy
	 *  as exp will yield inf.
	 */
	// if (r < 30) return (recursive_model (lamda, p, bc));
} // getThreshold

/* @brief	Exact recursive function to calculate Poisson Bionomial CDF

 Let S(n,r) be the probability of observing at a particular locus r errors
 in n reads at an error rate of p per base, p not necessarily identical.
 Then S can be defined by the following recurrence equations:
    Base case: S(0,0) = 1
    Induction: S(n,0) = S(n-1,0) * (1-p) when n>0
               S(n,n) = S(n-1,n-1) * p when n>0
               S(n,r) = S(n-1,r-1) * p + S(n-1,r) * (1-p) when 0<r<n
	log S(n,r): is log_nr in the code
 Dynamic programming can be calculated by storing only current column only
 0,0
 |   \
 1,0  1,1
 |   \ | \
 2,0 2,1  2,2

 We store log prob to prevent underflow problems using the fact that

    log(a+b) = log(1+exp(log(b)-log(a))) + log(a)

 For space and time efficiency, we calculate S(0,0) to S(n,0) for the base
 case, and then on each iteration i we use the values of the prev iter to
 calculate S(i,i) to S(n,i). Iterate until we find the lowest value i such
 that 1 - (S(n,0)+..+S(n,i))**bc < 0.05, where bc is the Sidak correction.
*/
int recursive_model (double lamda, const dvec_t& p, const double& bc) {
	int depth = p.size();

	std::vector<double> log_p (depth), log_q (depth), log_nr (depth);
	int negInf = -1000000;

	int r = 0;

	// calculate the first column starting from S(1,0)

    for (int i = 0; i < depth; ++ i) {
        log_p[i] = p[i] ? std::log(p[i]) : negInf;
        log_q[i] = (1 - p[i]) ? std::log(1 - p[i]) : negInf;
        log_nr[i] = !i ? log_q[i] : log_nr[i] = log_q[i] + log_nr[i-1];
    }

    double logsum = log_nr.back(); // log (S(n, 0))

    double totp = 1 - std::exp(logsum * bc);

    while (totp > 0.05 && r < depth) {

        double newlog_nr = log_p[r] + (r ? log_nr[r - 1] : 0); // log(s(n-1,r))

        for (int i = r + 1; i < depth; ++ i) {

            double newlog_q = newlog_nr + log_q[i];  //  log(s(n-1, r) * (1-p))
            double newlog_p = log_nr[i-1] + log_p[i]; // log(s(n-1,r-1) * p)
            double addedlog = (newlog_p > newlog_q) ? // log(s(n, r))
            		std::log(1 + std::exp(newlog_q - newlog_p)) + newlog_p :
            		std::log(1 + std::exp(newlog_p - newlog_q)) + newlog_q;

            log_nr[i-1] = newlog_nr;
            newlog_nr = addedlog;
        }

        logsum = std::log(1 + std::exp(newlog_nr - logsum)) + logsum;
        totp = 1 - std::exp(logsum * bc);
        ++ r;
    }

    {
    		std::cout << "recursive threshold = " << r + 1 << "\n";
    	}

    return (r + 1);
} // recursive_model

/* @brief Sidak correction	alpha <= 1 - (1-beta)^bc, where beta is the
 * 	individual test and alpha is the overal significance.
 * 	Derive the equation we have 1/bc (log (1 - alpha)) <= log (1 - beta)
 * 	Plug in CDF we have log (1 - alpha)/bc <= log F(X < x)
 *	Expand poisson we have
 *  log (1 - alpha)/bc + lamda  <= log \sum_{i=0}^{i=k-1} lamda^i / i!
 *	let [threshold] = log (1 - alpha)/bc + lamda
 *
 */
int poisson_model (const double& lamda, int bc, int n) {

	double threshold =  1.0/bc * std::log(0.95) + lamda;
	int k = 0;
	double sum = 1;
	double logsum = std::log(sum);
	/* direct calcualtion
	double numerator = 1; // lamda ^ k
	double factorial = 1; // k !
	while (logsum < threshold && k < n) {
		++ k;
		factorial *= k;
		numerator *= lamda;
		sum += (numerator/factorial);
		logsum = std::log(sum);
	}*/
	// log scale
	double factor = 0;
	while (logsum < threshold && k < n) {
		++ k;
		factor += std::log(k);
		sum += std::exp(k*std::log(lamda) - factor);
		logsum = std::log(sum);
	}

	return (k + 1);
} // poisson_model

/* @brief	Obtain Chi sqr p-value -- if certain entry of table is < 5,
 * 			using fisher's exact	test instead
 * Note: in [table] every row has identical size, so does the columns
 */
double chi_sqr_pvalue (const iivec_t& table){

	int num_row = table.size(), num_col;
	if (num_row < 2) return 1.0;
	num_col = table[0].size();
	if (num_col < 2) return 1.0;;

	int df = (num_row - 1) * (num_col - 1);
	// ----------- marginal table value ----------------
	ivec_t row_marginal (num_row, 0), col_marginal (num_col, 0);
	int sum = 0;
	for (int i = 0; i < num_row; ++ i) { // row
		for (int j = 0; j < num_col; ++ j) { // col
			row_marginal[i] += table[i][j];
			col_marginal[j] += table[i][j];
			sum += table[i][j];
		}
	}
	double chi_test_val = 0;
	bool is_test_proper = true;
	//--------------- expected table value -----------------
	ddvec_t exp_table (num_row, dvec_t (num_col));
	for (int i = 0; i < num_row; ++ i) { // row
		for (int j = 0; j < num_col; ++ j) { // col
			exp_table[i][j] = 1.0 * row_marginal[i] * col_marginal[j]/sum;
			chi_test_val += 1.0 * (exp_table[i][j] - table[i][j]) *
					(exp_table[i][j] - table[i][j])/exp_table[i][j];
			if (exp_table[i][j] < 5) is_test_proper = false;
		}
	}

	if (is_test_proper) {
		// this is how chi-squre p-value is implemented in boost
		return boost::math::gamma_q(df * 0.5, chi_test_val/ 2);
	} else return -1.0;
}

/* @brief
 *
 *			x00	 x01		a
 *			x10	 x11		b
 *			c	 d      sum
 */
double tbt_fisher_exact_twotail_pval (int x00, int x01, int x10, int x11){
	static LogFac lf;

	int a = x00 + x01, b = x10 + x11, c = x00 + x10, d = x01 + x11;
	int sum = c + d;
	int min_x00 = x00 - std::min(x00, x11),
		max_x00 = std::min (a, c);
	dvec_t all_p;
	double critical_value = x00;
	for (x00 = min_x00; x00 <= max_x00; ++ x00) {
		x01 = a - x00;
		x10 = c - x00;
		x11 = d - x01;
		double log_geo = lf.get_value(a) + lf.get_value(b) +
			lf.get_value(c) + lf.get_value(d) - lf.get_value(x00) -
			lf.get_value(x01) - lf.get_value(x10) - lf.get_value(x11) -
			lf.get_value(sum);
		all_p.push_back(std::exp(log_geo));
		if (critical_value == x00) critical_value = all_p.back();
		//sum_log_pval += std::exp(log_geo);
		//std::cout << std::exp(log_geo) << "\n";
	}
	double pval = 0.0;
	for (int i = 0; i <= (int) all_p.size(); ++ i) {
		if (all_p[i] <= critical_value) pval += all_p[i];
	}
	return pval;
}

/* @brief	CDF \phi(X) of Normal distribution; code taken rom
 * 			http://www.johndcook.com/cpp_phi.html
 */
/*
double cdf_normal(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
} // cdf_normal
*/
/* @brief	Srinivasa Ramanujan ln(n!) factorial estimator
 */
/*
double factorial_approx (int n) {
	if (n < 2) return 0;
	double a = n * log(n) - n;
	double b = log(n * (1 + 4 * n * (1 + 2 * n))) / 6;
	return a + b + log(3.14159265358979323846) / 2;
}
*/
/* @brief CDF of Poisson distribution: F(x <= k)
 */
/*
double cdf_poisson (double lamda, int k) {
	double sum_log_factorial_k[21] = {0, 0, 0.693147181, 1.79175947,
		3.17805383, 4.78749174,  6.57925121, 8.52516136, 10.6046029,
		12.8018275,  15.1044126, 17.5023078,  19.9872145, 22.5521639,
		25.1912212, 27.8992714, 30.6718601, 33.5050735, 36.3954452,
		39.3398842, 42.3356165};
	// e^{-lamda}
	double cdf = exp(-1*lamda);
	// 1 + \sum_{i=1}^{i=k} exp {iln(lamda) - (lni!)}
	double sum = 1.0;
	for (int i = 1; i <= k; ++ i) {
		double logs = 0.0;
		if (i > 20) logs = factorial_approx(i);
		else logs = sum_log_factorial_k[i];
		sum += std::exp(i*log(lamda) - logs);
	}
	cdf *= sum;
	return cdf;
} // cdf_poisson
*/
