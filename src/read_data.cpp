/*
 * read_data.cpp
 *
 *  Created on: Sep 19, 2011
 *      Author: bkloppenborg
 */

#include "read_data.h"
#include <cstdlib>
#include <stdexcept>

#include "ReadTextFile.h"

void read_data_ast(string filename, string comment_chars, vector< vector<int> > split_info, vector<vector<double> > data, bool read_no_error, double default_error)
{
	// First determine the type of observation and fork it off to the
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Astrometry Data File");
	vector < vector<double> > results;

	double t, x, y, e_x, e_y, P_alpha, P_delta;
	double r, theta;
	int i_x = 1;
	int i_ex = 2;
	int i_y = 3;
	int i_ey = 4;
	int i_Pa = 5;
	int i_Pd = 6;

	if(read_no_error)
	{
		i_y -= 1;
		i_ex = 5;
		i_ey = 5;
		i_Pa -= 2;
		i_Pd -= 2;

	}

	// Now iterate through the lines and tokenize them.  Notice, we use vector.at(n) instead of vector[n]
	// so that we can catch the signals from exceptions.  The data file is only parsed once so this is ok.

	for (unsigned int i = 0; i < lines.size(); i++)
	{
		vector < string > tokens;

		if(split_info.size() > 0)
			tokens = Tokenize(lines[i], split_info);
		else
			tokens = Tokenize(lines[i]);

		// And now attempt to read in the line
		try
		{
			t = atof(tokens.at(0).c_str());
			x = atof(tokens.at(i_x).c_str());
			y = atof(tokens.at(i_y).c_str());
		}
		catch(...)
		{
			throw std::runtime_error("Could not parse line in astrometric data file.");
		}

		// Now for the uncertainties (these may not exist)
		try
		{
			e_x = atof(tokens.at(i_ex).c_str());
		}
		catch(...)
		{
			e_x = 0;
		}

		try
		{
			e_y = atof(tokens.at(i_ey).c_str());
		}
		catch(...)
		{
			e_y = 0;
		}

		// Lastly the parallax factors
		try
		{
			P_alpha = atof(tokens.at(i_Pa).c_str());
		}
		catch(...)
		{
			P_alpha = 0;
		}

		try
		{
			P_delta = atof(tokens.at(i_Pd).c_str());
		}
		catch(...)
		{
			P_delta = 0;
		}

		if(e_x == 0)
			e_x = default_error;
		if(e_y == 0)
			e_y = default_error;

		// Enable if you want to see the data.
		//printf("%i %i %i %i %i %i %i \n", 0, i_x, i_ex, i_y, i_ey, i_Pa, i_Pd);
		//printf("%f %f %f %f %f %f %f \n", t, x, e_x, y, e_y, P_alpha, P_delta);

		// Push this station on to the list of stations for this array.
		vector<double> temp;
		temp.push_back(t);
		temp.push_back(x);
		temp.push_back(e_x);
		temp.push_back(y);
		temp.push_back(e_y);
		temp.push_back(P_alpha);
		temp.push_back(P_delta);
		data.push_back(temp);
	}
}
