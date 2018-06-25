// Copyright (c) 2012 Andre Martins
// All Rights Reserved.
//
// This file is part of AD3 2.1.
//
// AD3 2.1 is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// AD3 2.1 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with AD3 2.1.  If not, see <http://www.gnu.org/licenses/>.

#ifndef UTILS_H_
#define UTILS_H_

#if defined(_WIN32) && !defined(__MINGW32__)
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <vector>
#include <string>
#include <algorithm>

#if defined(_WIN32) && !defined(__MINGW32__)
#include <windows.h> //I've ommited this line.
#ifndef _WINSOCKAPI_
struct timeval {
        long    tv_sec;         /* seconds */
        long    tv_usec;        /* and microseconds */
};
#endif
extern int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif

namespace AD3 {

#define NEARLY_EQ_TOL(a,b,tol) (((a)-(b))*((a)-(b))<=(tol))
#define NEARLY_BINARY(a,tol) (NEARLY_EQ_TOL((a),1.0,(tol)) || NEARLY_EQ_TOL((a),0.0,(tol)))
#define NEARLY_ZERO_TOL(a,tol) (((a)<=(tol)) && ((a)>=(-(tol))))

#define MAX(a,b) (((a)<(b))? (b) : (a))

#define LOG_ZERO -1000
#define LOG_STABLE(a) (a == 0? LOG_ZERO : log(a))

extern int diff_ms(timeval t1, timeval t2);

extern int diff_us(timeval t1, timeval t2);

extern void InsertionSort(std::pair<double, int> arr[], int length);

extern int project_onto_simplex_cached(double* x,
				       int d,
				       double r, 
				       std::vector<std::pair<double,int> >& y);

extern int project_onto_simplex(double* x, int d, double r);

extern int project_onto_cone_cached(double* x, int d,
				    std::vector<std::pair<double,int> >& y);
				    
extern int project_onto_budget_constraint(double* x, int d, double budget);	

extern int project_onto_budget_constraint_cached(double* x,
                                                 int d,
                                                 double budget, 
                                                 std::vector<std::pair<double,int> >& y);

extern int project_onto_knapsack_constraint(double* x, double* costs, int d,
                                            double budget);

extern int solve_canonical_qp_knapsack(const std::vector<double> &lower_bounds,
                                       const std::vector<double> &upper_bounds,
                                       const std::vector<double> &weights,
                                       double total_weight,
                                       std::vector<double> *solution);

                                                 
extern void StringSplit(const std::string &str,
			const std::string &delim,
			std::vector<std::string> *results);

extern void TrimComments(const std::string &delim, std::string *line);

extern void TrimLeft(const std::string &delim, std::string *line);

extern void TrimRight(const std::string &delim, std::string *line);

extern void Trim(const std::string &delim, std::string *line);

} // namespace AD3

#endif
