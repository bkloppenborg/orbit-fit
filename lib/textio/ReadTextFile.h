/*
 * ReadTextFile.h
 *
 *  Created on: Feb 18, 2011
 *      Author: bkloppenborg
 */

#ifndef READTEXTFILE_H_
#define READTEXTFILE_H_

#include <vector>
#include <string>
using namespace std;

void 	StringSplit(string str, string delim, vector<string> & results);
string 	StripWhitespace(string str);
void 	StripWhitespace(vector<string> & strings);
vector<string> ReadFile(string filename, string comment_chars, string error_message);

vector<string> Tokenize(string line);
vector<string> Tokenize(string line, vector< vector<int> > split_locations);

bool 	FileExists(string filename);

#endif /* READTEXTFILE_H_ */
