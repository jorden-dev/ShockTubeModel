#include <stdio.h>
#include <math.h>
#include "MatrixOps.h"

// LU Decompose with partial pivoting
//      Inputs: A - Matrix of dimension NxN, Tol - detect
//          faulure of matrix if degenerate, P - Permutation
//          array of size N+1.
//      Outputs: A is changed to contain both L and U such
//          that P*A=L*U.
int LUDecompose(int N, double A[N][N], double Tol, int *P)
{
    int i, j, k, imax;
    double maxA, ptr, absA;

    // Unit permutation matrix P[N]
    for (int i = 0; i <= N; i++)
    {
        P[i] = i;
    }

    for (int i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        // Find row with max value
        for (int k = i; k < N; k++)
        {
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }
        }

        // Check matrix
        if (maxA < Tol)
        {
            printf("%s","Failure, matrix is degenerate");
            return 0;
        }

        // Pivot matrices
        if (imax != i)
        {
            // pivot P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            // pivot rows of A
            // TODO: pivot using pointers
            for (j = 0; j < N; j++)
            {
                ptr = A[i][j];
                A[i][j] = A[imax][j];
                A[imax][j] = ptr;
            }
            
            // count pivots starting from N
            P[N]++;
        }

        for (j = i + 1; j < N; j++)
        {
            A[j][i] /= A[i][i];
            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
        
    }
   return 0;
}

// Determinant
// TODO: Debug!
double determinant(int N, double A[N][N], int *P)
{
    double det = A[0][0];

    for (int i = 1; i < N; i++)
    {
        det *= A[i][i];
    }
    return (P[N] - N);
}

// Invert 
double invert(int N, double A[N][N], int *P, double IA[3][3])
{
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            IA[i][j] = P[i] == j;
            for(int k = 0; k < i; k++)
            {
                IA[i][j] -= A[i][k] * IA[k][j];
            }
        }

        for (int i = N -1; i >= 0; i--)
        {
            for (int k = i + 1; k < N; k++)
            {
                IA[i][j] -= A[i][k] * IA[k][j];
            }
            IA[i][j] /= A[i][i];
        }
    }
}

// Multiply two matrices
void multiply(int r1, int c1, int r2, int c2,
                double A[r1][c1],
                double B[r2][c2],
                double C[r1][c2])
{
    // Check matrices A and B for correct dimensions
    if (c1 != r2)
    {
        printf("Invalid matrice dimensions for multiplication!");
        return;
    }

    // Set C elements to zero
    for (int i = 0; i < r1; ++i)
    {
        for (int j = 0; j < c2; ++j)
        {
            C[i][j] = 0;
        }
    }
    
    for (int i = 0; i < r1; ++i)
    {
        for (int j = 0; j < c2; ++j)
        {
            for (int k = 0; k < c1; ++k)
            {
                C[i][j] += A[i][k] * B[k][j];
            } 
        }
    }
}

// Multiply matrix with vector
// Inputs: A = r1xc1 matrix, B = r2x1 matrix
void multiplyV(int r1, int c1, int r2, int c2,
                double A[r1][c1],
                double B[r2],
                double C[r1][c2])
{
    // Check matrices A and B for correct dimensions
    if (c1 != r2)
    {
        printf("Invalid matrice dimensions for multiplication!");
        return;
    }
    
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            for (int k = 0; k < c1; k++)
            {
                C[i][j] += A[i][k] * B[k];
            } 
        }
    }
}

// Set matrix elements to zero
void zeroMatrix(int N, double A[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = 0;
        }
        
    }
    
}

// Transpose

// Linear spaced vector
//      Inputs: A - Matrix of dimension N, start - first element
//          value, end - final element value.
//      Outputs: A is changed to be a linearspaced vector 
//          start:N:end.
void linspace(int N, double A[N], double start, double end)
{
    for (int i = 0; i < N; i++)
    {
        A[i] = i/(N - end);
    }
}

// Max vector element
double vectormax(int N, double A[N]){
    int index = 0;
    for (int i = 0; i < N; i++){
        if (A[i] > A[index]){
            index = i;
        }
    }
    return index;
}
