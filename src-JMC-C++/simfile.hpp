/*
 * simfile.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#ifndef SIMFILE_H_
#define SIMFILE_H_

void dosimfile(int seed);
void dosimstat(int seed);
int detpstarOK(int nsel, int scen, long double **phistar);

#endif /* SIMFILE_H_ */
