/*
 * fitboth.h
 *
 *  Created on: Sep 6, 2011
 *      Author: bkloppenborg
 */

#ifndef FITBOTH_H_
#define FITBOTH_H_

void print_help();

void run_fit();

void ParseProgOptions(int argc, char *argv[], bool & param_error);

int main(int argc, char *argv[]);

#endif /* FITBOTH_H_ */
