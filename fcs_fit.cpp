/* fcs-tools - Tools for FCS data analysis
 *
 * Copyright © 2010 Ben Gamari
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 *
 * Author: Ben Gamari <bgamari@physics.umass.edu>
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <utility>
#include <cassert>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_multifit_nlin.h>


using std::vector;
using std::string;
using std::array;

struct point {
	uint64_t time;
	double value;
};

/*
 * data_set: Represents a dataset and its corresponding physical parameters
 *
 * While really only num_density will vary over experiments (in most cases), 
 * special cases are evil. Thus we make  all physical parameters per-dataset,
 * despite diff_time and aspect_ratio being constant across sets.
 */
struct data_set {
	// Number density (N)
	double num_density;
	// Characteristic diffusion time (tau)
	double diff_time;
	// Measurement volume aspect ratio (R)
	double aspect_ratio;
	// Measured correlation function
	vector<point> points;
};


/*
 * diffusion_fit_f: Evaluate G(tau)
 */
double diffusion_fit_f(const gsl_vector* x, vector<point>& points,
		double num_density, double diff_time, double aspect_ratio) {
	double tau_taud = gsl_vector_get(x,0) / diff_time;
	double b = (1 + tau_taud) / (1 + pow(aspect_ratio, -2) * tau_taud);
	return (1 / num_density / sqrt(b));
}

/*
 * diffusion_fit_df: Evaluate (dG/dN, dG/dtau, dG/dR)
 */
array<double,3> diffusion_fit_df(const gsl_vector* x, vector<point>& points,
		double num_density, double diff_time, double aspect_ratio) {
	array<double,3> res;
	// dG/dN
	res[0] = 0;
	// dG/dtau
	res[1] = 0;
	// dG/dR
	res[2] = 0;
}


int fit_func_f(const gsl_vector* x, void* _data, gsl_vector* f) {
	vector<data_set*>* data = (vector<data_set*>*) _data;
	for (int i=0; i<data->size(); i++) {
		data_set* ds = (*data)[i];
		double v = diffusion_fit_f(x, ds->points,
				ds->num_density, ds->diff_time, ds->aspect_ratio);
		gsl_vector_set(f, i, v);
	}
	return GSL_SUCCESS;
}

int fit_func_df(const gsl_vector* x, void* _data, gsl_matrix* J) {
	vector<data_set*>* data = (vector<data_set*>*) _data;
	for (int i=0; i<data->size(); i++) {
		data_set* ds = (*data)[i];
		array<double,3> df = diffusion_fit_df(x, ds->points,
				ds->num_density, ds->diff_time, ds->aspect_ratio);
		for (int j=0; j<3; j++)
			gsl_matrix_set(J, j, i, df[j]);
	}
	return GSL_SUCCESS;
}

int fit_func_fdf(const gsl_vector* x, void* _data, gsl_vector* f, gsl_matrix* J) {
	vector<data_set*>* data = (vector<data_set*>*) _data;
	for (int i=0; i<data->size(); i++) {
		data_set* ds = (*data)[i];
		double v = diffusion_fit_f(x, (*data)[i]->points,
				ds->num_density, ds->diff_time, ds->aspect_ratio);
		gsl_vector_set(f, i, v);

		array<double,3> df = diffusion_fit_df(x, ds->points,
				ds->num_density, ds->diff_time, ds->aspect_ratio);
		for (int j=0; j<3; j++)
			gsl_matrix_set(J, j, i, df[j]);

	}
	return GSL_SUCCESS;
}

void fit(vector<data_set*>& data) {
	const size_t n_params = 3;
	const size_t n_sets = data.size();
	const double epsabs = 1e-4, epsrel = 1e-4;

	int res;
	gsl_multifit_fdfsolver* solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder,
			data.size(), n_params);
	double x_init[n_params] = { 0, };
	gsl_vector_view x = gsl_vector_view_array(x_init, n_params);

	gsl_multifit_function_fdf f;
	f.f = &fit_func_f;
	f.df = &fit_func_df;
	f.fdf = &fit_func_fdf;
	f.n = n_sets;
	f.p = n_params;
	f.params = &data;

	res = gsl_multifit_fdfsolver_set(solver, &f, &x.vector);
	assert(res == GSL_SUCCESS);

	bool converged = false;
	while (!converged) {
		gsl_multifit_fdfsolver_iterate(solver);
		res = gsl_multifit_test_delta(solver->dx, solver->x, epsabs, epsrel);
		converged = (res == GSL_SUCCESS);
	}

	gsl_multifit_fdfsolver_free(solver);
}

int main(int argc, char** argv) {
	vector<data_set*> data;

	string line;
	while (std::getline(std::cin, line, '\n')) {
		data_set* ds = new data_set();
		boost::tokenizer<> tokens(line);
		auto tok = tokens.begin();
		string file = *tok; tok++;
		ds->num_density = boost::lexical_cast<double>(*tok); tok++;
		ds->diff_time = boost::lexical_cast<double>(*tok); tok++;
		ds->aspect_ratio = boost::lexical_cast<double>(*tok);

		std::ifstream f(file);
		while (!f.eof() && !f.fail()) {
			point p;
			std::cin.read((char*) &p.time, sizeof(uint64_t));
			std::cin.read((char*) &p.value, sizeof(double));
			ds->points.push_back(p);
		}
	} while (!std::cin.fail() && !std::cin.eof());

	fit(data);
}

