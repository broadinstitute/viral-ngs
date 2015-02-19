//========================================================================
// Project     : VariantCaller
// Name        : stat_func.h
// Author      : Xiao Yang
// Created on  : Sep 25, 2012
// Version     : 1.0
// Copyright Broad Institute, Inc. 2013.
// Notice of attribution: The V-Phaser 2.0 program was made available through the generosity of Genome Sequencing and Analysis Program at the Broad Institute, Inc. per Yang X, Charlebois P, Macalalad A, Henn MR and Zody MC (2013) V-Phaser 2.0: Variant Inference for Viral Populations” See accompanying file LICENSE_1_0.txt.  Distribution subject to licenses from Boost Software and MIT (http://www.boost.org/LICENSE_1_0.txt and https://github.com/pezmaster31/bamtools/blob/master/LICENSE).

// Description :
//========================================================================


#ifndef STAT_FUNC_H_
#define STAT_FUNC_H_

#include "xutil.h"
//#include "boost/math/distributions/chi_squared.hpp"
#include <boost/math/distributions/gamma.hpp>
/* stores vector of log factorial
 */
class LogFac {
private:
	dvec_t	_log_prefix_sum;
public:
	LogFac () {
		_log_prefix_sum.push_back(0); // initialize log 0 !
		_log_prefix_sum.push_back(0); // initialize log 0 !
	};
	double get_value (int n) { // log (n!)
		int sz = _log_prefix_sum.size();
		for (int i = sz; i <= n; ++ i) {
			_log_prefix_sum.push_back(_log_prefix_sum[i - 1] + std::log(i));
		}
		return _log_prefix_sum[n];
	}
};

int getThreshold (double lamda, const dvec_t& p, int bc);

int poisson_model(const double& lamda, int bc, int n);
int recursive_model (double lamda, const dvec_t& p, const double& bc);

double chi_sqr_pvalue (const iivec_t& table);

double tbt_fisher_exact_twotail_pval (int x00, int x01, int x10, int x11);

double cdf_poisson (double lamda, int k);

#endif /* STAT_FUNC_H_ */
