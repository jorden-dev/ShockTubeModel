#ifndef MATRIXOPS_H_
#define MATRIXOPS_H _
// Various functions to work with matrix operations.

// LU Decompose
int LUDecompose(int N, double A[N][N], double Tol, int *P);

// Determinant
double determinant(int N, double A[N][N], int *P);

// Invert 
double invert(int N, double A[N][N], int *P, double IA[3][3]);

// Transpose

// Multiply
void multiply(int r1, int c1, int r2, int c2,
                double A[r1][c1],
                double B[r2][c2],
                double C[r1][c2]);

// Multiple nxn matrix by 1xn vector
void multiplyV(int r1, int c1, int r2, int c2,
                double A[r1][c1],
                double B[r2],
                double C[r1][c2]);

// Set matrix elements to zero
void zeroMatrix(int N, double A[N][N]);

// Linear spaced vector
void linspace(int N, double A[N], double start, double end);

// Max vector element
double vectormax(int N, double A[N]);

#endif