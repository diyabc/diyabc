// matrices.h
//

#ifndef MATRICES_H_
#define MATRICES_H_

#define EPSILON 9.62964972193617926528E-35

void libereD (int n, double * * A);
void libereL (int n, long double * * A);
double * * transpose (int n, int m, double * * A);
double * * prodM (int n, int m, int p, double * * A, double * * B);
double * * prodMs (int n, int m, double * * A, double b);
long double * * transposeL (int n, int m, long double * * A);
long double * * prodML (int n, int m, int p, long double * * A, long double * * B);
long double * * prodMsL (int n, int m, long double * * A, long double b);
double * * invM (int n, double * * A);
int inverse (int n, long double * * A, long double * * C);
int jacobi (int n, double * * A, double * D, double * * V);
int jacobiL (int n, long double * * A, long double * D, long double * * V);
double kappa (int n, long double * * A);
int inverse_Tik (int n, long double * * A, long double * * C);
int inverse_Tik2 (int n, long double * * A, long double * * C, long double coeff);

#endif
