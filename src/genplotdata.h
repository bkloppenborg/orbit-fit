/*
 * genplotdata.h
 *
 *  Created on: Sep 19, 2011
 *      Author: bkloppenborg
 */

#ifndef PARSE_RESULTS_H_
#define PARSE_RESULTS_H_

#include <vector>
#include <string>
using namespace std;

// Now for the functions

void extract_params_ast(vector<double> params,
		double & omega, double & e, double & tau, double & T,
		double & Omega, double & inc, double & alpha,
		double & x0, double & y0, double & mu_x, double & mu_y, double & pi);

void extract_params_both(vector<double> & params, bool & fit_pm,
		double & omega, double & e, double & tau, double & T,
		double & K, double & gamma,
		double & Omega, double & inc, double & alpha,
		double & x0, double & y0, double & mu_x, double & mu_y, double & pi);

void extract_params_rv(vector<double> params, bool & fit_pm,
		double & omega, double & e, double & tau, double & T,
		double & K, double & gamma);

int main(int argc, char *argv[]);

void print_help();

void gendata_ast(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step);

void gendata_ast(string output_name, bool fit_pm,
		double omega, double e, double tau, double T, double Omega, double inc, double alpha,
		double x0, double y0, double mu_x, double mu_y, double pi,
		double t_min, double t_max, double t_step);

void gendata_both(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step);

void gendata_rv(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step);

void gendata_rv(string output_name,
		double omega, double e, double tau, double T,
		double K, double gamma,
		double t_min, double t_max, double t_step);

void residuals_ast(string output_basename, vector<vector<double> > & params, vector<vector<double> > & data);

void residuals_ast(string output_name, bool fit_pm,
		double omega, double e, double tau, double T,
		double Omega, double inc, double alpha,
		double x0, double y0, double mu_x, double mu_y, double pi,
		vector<vector<double> > & data);

void residuals_both(string output_basename, vector<vector<double> > & params, vector<vector<double> > & ast_data, vector<vector<double> > & rv_data);

void residuals_rv(string output_name,
		double omega, double e, double tau, double T,
		double K, double gamma,
		vector<vector<double> > & data);

void residuals_rv(string output_basename, vector<vector<double> > & params, vector<vector<double> > & data);

vector<vector<double> > read_data_summary(string filename, string comment_chars);

#endif /* PARSE_RESULTS_H_ */
