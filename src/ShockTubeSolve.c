/*
Program to simultae non-dimensional 1D shock tube using different Euler schemes.

Author: Jorden Schulte 
Assignment #3

*/

#include <stdio.h>
#include "MatrixOps.h"
#include "Solve.h"
#include <math.h>
#include <string.h>


int main(){
    /* ------------------------------- USER INPUTS ------------------------------ */
    // "_L" designates left of diaphram, "_R" is right of diaphram
    double p_L = 1e5;               // Pressure [Pa]
    double rho_L = 1;               // Density [kg/m3]
    double u_L = 0;                 // Velocity [m/s]

    double p_R = 1e3;
    double rho_R = 0.01;
    double u_R = 0;

    int nx = 101;                 // Number of grid points 
    double x0 = 0.5;                // Inital diaphragm position
    double xmin = 0;                // Front of plate
    double xmax = 1.0;              // End of plate

    double cfl = 0.9;               // CFL condition
    double runtime = 4.9e-4;        // Solver runtime
    double starttime = 0.0;         // Solver start time
    // NOTE: For production code inputs would be read from input file.

    /* ----------------------------- INITIALIZE CASE ---------------------------- */
    // Construct grid
    double dx = (xmax - xmin) / nx;
    double x[nx];
    linspace(nx, x, xmin, xmax);

    // Initialize U and F(U) matrices of conservative vectors
    double U[3][nx];
    double FU[3][nx];
    double p[nx];
    double u[nx];
    double E[nx];
    double c[nx];

    double E_L, E_R; // Constant Energy in L and R sections
    E_L = perfgas_e(p_L, rho_L, u_L);
    E_R = perfgas_e(p_R, rho_R, u_R);

    double c_L, c_R; // Speed of sound
    c_L = sqrt(GAMMA * p_L / rho_L);
    c_R = sqrt(GAMMA * p_R / rho_R);

    for (int i = 0; i < nx; i++){
        if (x[i] < x0){
            U[0][i] = rho_L;
            U[1][i] = rho_L * u_L;
            U[2][i] = rho_L * E_L;

            FU[0][i] = rho_L * u_L;
            FU[1][i] = rho_L * pow(u_L, 2) + p_L;
            FU[2][i] = rho_L * u_L * (E_L + p_L / rho_L);

            u[i] = u_L;
            c[i] = c_L;
            p[i] = p_L;
        }
        if (x[i] >= x0){
            U[0][i] = rho_R;
            U[1][i] = rho_R * u_R;
            U[2][i] = rho_R * E_R;

            FU[0][i] = rho_R * u_R;
            FU[1][i] = rho_R * pow(u_R, 2) + p_R;
            FU[2][i] = rho_R * u_R * (E_R + p_R / rho_R);

            u[i] = u_R;
            c[i] = c_R;
            p[i] = p_R;
        }
    }

    // Max wave speed
    int maxindex;
    double ucmax;
    double uc[nx];

    /* ----------------------------------- RUN ---------------------------------- */
    double timenow;         // Current solver iteration time
    int iter = 1;           // Iteration counter
    double dt;    // Delta t, timestep

    // Begin time step loop
    while (timenow < runtime){

        // Calculate stable dt from cfl
        for (int i = 0; i<nx; i++){
            uc[i] = u[i] + c[i];
        }
        maxindex = vectormax(nx, uc);
        ucmax = uc[maxindex];
        dt = dtcalc(cfl, dx, ucmax);

        // Solve at current time step
        // Update time
        timenow += dt;
        iter++;

        // Solver scheme choices:
        // "centraldiff" = Central difference with 1st order scalar dissipation.
        // "laxwendroff" = CD in conservative Lax-Wendroff form.
        // "vanleer" = Van Leer flux-vector splitting.
        // "maccormick" = MacCormick corrector-predictor scheme.

        //centraldiff(nx, U, FU, u, p, c, dt, dx);
        laxwendroff(nx, U, FU, u, p, c, dt, dx);
        //vanleer(nx, U, FU, u, p, c, dt, dx);
        //maccormack(nx, U, FU, u, p, c, dt, dx);

        // Terminal Output
        printf("\nTime = %9.6f\tIter. = %d", timenow, iter);

    }

    // Postprocess other quantities for plotting:
    // Density, Mach number, mass flow, entropy (approx).
    double rho[nx], M[nx], m[nx], s[nx];
    for (int i = 0; i < nx; i++){
        rho[i] = U[0][i];
        M[i] = u[i]/c[i];
        m[i] = rho[i] * u[i];
        s[i] = p[i] / pow(rho[i], GAMMA) - p_L / pow(rho_L, GAMMA);
        s[i] = s[i]/RAIR;
    }

    /* ----------------------------- OUTPUT RESULTS ----------------------------- */
    char filename[100] = "ShockTubeSolve_VL"; // TODO: Use scheme choice for file name
    strcat(filename, "_101"); // TODO: Convert int nx to string for file name
    strcat(filename, ".csv");
    printf("\nCreating %s file", filename);

    FILE *fp;
    fp = fopen(filename, "w+");

    // Column headers
    int n_cols = 9;
    // NOTE: y and z coord. not needed here, added for paraview compatability
    double y[nx], z[nx];
    fprintf(fp, "x coord, y coord, z coord, p, s/R, u, M, rho, m");
    // Column values
    for (int i = 0; i < nx; i++)
    {
        fprintf(fp, "\n");

        fprintf(fp, "%9.6f,\t", x[i]);
        fprintf(fp, "%9.6f,\t", y[i]);
        fprintf(fp, "%9.6f,\t", z[i]);
        fprintf(fp, "%9.6f,\t", p[i]);
        fprintf(fp, "%9.6f,\t", s[i]);
        fprintf(fp, "%9.6f,\t", u[i]);
        fprintf(fp, "%9.6f,\t", M[i]);
        fprintf(fp, "%9.6f,\t", rho[i]);
        fprintf(fp, "%9.6f,\t", m[i]);

    }
    fclose(fp);
    printf("\nFile %s created\n", filename);
    return 0;
}