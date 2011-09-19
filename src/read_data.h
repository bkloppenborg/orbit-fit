/*
 * read_data.h
 *
 *  Created on: Sep 19, 2011
 *      Author: bkloppenborg
 */

#ifndef READ_DATA_H_
#define READ_DATA_H_

#include <vector>
#include <string>

using namespace std;

void read_data_ast(string filename, string comment_chars, vector< vector<int> > split_info, vector<vector<double> > & data, bool read_no_error, double default_error);

void read_data_rv(string filename, string comment_chars, vector< vector<int> > split_info, vector<vector<double> > & data, bool read_no_error, double default_error);

#endif /* READ_DATA_H_ */
