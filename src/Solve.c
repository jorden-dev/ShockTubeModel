#include "Solve.h"
#include <math.h>
#include <string.h>
#include "MatrixOps.h"

const double GAMMA = 1.4;
const double RAIR = 287.05;

/* ------------------------------- Perfect Gas ------------------------------ */
// Solves for pressure using perfect gas relations.
// Inputs: E = total energy, rho = density, u = velocity scalar.
// Outputs: p = pressure.
double perfgas_p(double E, double rho, double u){
    double p;
    p = (GAMMA - 1) * (rho * E - 0.5 * rho * pow(u, 2));
    return p;
}

// Solves for energy using perfect gas relation.
// Inputs: p = pressure, rho = density, u = scalar velocity.
double perfgas_e(double p, double rho, double u){
    double E;
    E = p / (rho * (GAMMA - 1)) + 0.5 * pow(u, 2);
    return E; 
}

/* ------------------------------ CFL Condition ----------------------------- */
// INPUTS: cfl = target cfl number, dx = grid spacing (assume uniform), umax = 
// max wave speed.
// OUTPUTS: dt = Delta t timestep.
double dtcalc(double cfl, double dx, double umax){
    return dx * cfl / umax;
}

/* ----------------------------- Flux Functions ----------------------------- */
// Flux per Lax-Wendroff
void flux_lw(int N, double flux[3][N-1], double U[3][N], double FU[3][N], double u[N], double c[N], double residual[3][N]){
    // Evaluate u + c at face
    double ucface[N-1];
    for (int i = 0; i < N-1; i++){
        //ucface[i] = fabs(0.5 * (u[i+1] + u[i])) + 0.5 * (c[i+1] + c[i]);
        ucface[i] = 0.5 * ((fabs(u[i]) + c[i]) + (fabs(u[i+1]) + c[i+1]));
    }

    // Loop over faces
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            flux[j][i] = 0.5 * (FU[j][i+1] + FU[j][i]) 
                        - 0.5 * ucface[i] * (U[j][i+1] - U[j][i]);
            residual[j][i] = residual[j][i] + flux[j][i];
            residual[j][i+1] = residual[j][i+1] - flux[j][i];
        }
    } 
}

// Flux vector splitting per Vanleer
void flux_vl(int N, double flux[3][N-1], double U[3][N], double FU[3][N], double u[N], double c[N], double residual[3][N]){
    // Evaluate Mach in each cell
    double M[N];
    double fplus[3][N];
    double fminus[3][N];
    for (int i = 0; i < N; i++){
        M[i] = u[i] / c[i];

        if (M[i] > 1){
            fplus[0][i] = FU[0][i];
            fplus[1][i] = FU[1][i];
            fplus[2][i] = FU[2][i];
            fminus[0][i] = 0;
            fminus[1][i] = 0;
            fminus[2][i] = 0;
        }
        else if (M[i] <= 1 && M[i] >= -1){
            fplus[0][i] = 0.25*U[0][i]*c[i] * pow((M[i] + 1),2) * 1;
            fplus[1][i] = 0.25*U[0][i]*c[i] * pow((M[i] + 1),2)
                        * 2 * c[i] / GAMMA * (1 + 0.5*(GAMMA-1)*M[i]);
            fplus[2][i] = 0.25*U[0][i]*c[i] * pow((M[i] + 1),2)
                        * 2*pow(c[i],2) / (pow(GAMMA,2) - 1) 
                        * pow((1 + 0.5*(GAMMA-1)*M[i]),2);
            
            fminus[0][i] = -0.25*U[0][i]*c[i] * pow((M[i] - 1),2) * 1;
            fminus[1][i] = -0.25*U[0][i]*c[i] * pow((M[i] - 1),2) 
                        * 2 * c[i] / GAMMA * (-1 + 0.5*(GAMMA-1)*M[i]);
            fminus[2][i] = -0.25*U[0][i]*c[i] * pow((M[i] - 1), 2)
                        * 2*pow(c[i], 2) / (pow(GAMMA,2) - 1) 
                        * pow((1 - 0.5*(GAMMA-1) * M[i]), 2);
        }
        else if(M[i] < -1){
            fplus[0][i] = 0;
            fplus[1][i] = 0;
            fplus[2][i] = 0;
            fminus[0][i] = FU[0][i];
            fminus[1][i] = FU[1][i];
            fminus[2][i] = FU[2][i];   
        }
    }

    // Determine flux through faces f[i - 1/2]
    // Loop over faces
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            flux[j][i] = fplus[j][i] + fminus[j][i+1];
            residual[j][i] = residual[j][i] + flux[j][i];
            residual[j][i+1] = residual[j][i+1] - flux[j][i];
        }
    }
}

/* ---------------------------- Update Flow Field --------------------------- */
void updateflow(int N, double Unew[3][N], double U[3][N], double FU[3][N], 
                double u[N], double c[N], double p[N]){
    // Update U, u, p, and c
    for (int i = 1; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            U[j][i] = Unew[j][i];
        }

        u[i] = U[1][i] / U[0][i];
        p[i] = perfgas_p(U[2][i]/U[0][i], U[0][i], u[i]); 
        c[i] = sqrt(GAMMA * p[i] / U[0][i]);
    }

    // Update F(U)
    for (int i =1; i < N-1; i++){
        FU[0][i] = U[1][i];
        FU[1][i] = U[1][i] * u[i] + p[i];
        FU[2][i] = U[1][i] * (U[2][i] + p[i])/U[0][i];
    }
}

/* --------------------------- Boundary Conditions -------------------------- */

/* --------------------------- Central Difference --------------------------- */
// Inputs: N = number of grid points, U = vector of conserved variables where
// U[0][:] = rho, U[1][:] = rho*u, U[2][:] = rho*E, FU = vector of flux vectors
// where FU[0][:] = rho*u, F[1][:] = rho*u^2+p, FU[2][:] = rho*u*H,
// u = primitive velocity, c = speed of sound, dt = delta time, dx = delta x.
// Outputs: Updates U, FU, u, p, and c for current timestep.

void centraldiff(int N, double U[3][N], double FU[3][N], double u[N], double p[N],
                double c[N], double dt, double dx){ 
    // Initialize temporary U_new matrix
    double Unew[3][N];

    // Loop over points to update
    for (int i = 1; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            Unew[j][i] = U[j][i] - dt/dx * (0.5*(FU[j][i+1] - FU[j][i-1])
                        - 0.5*(fabs(u[i]) + c[i]) 
                        * (U[j][i+1] - 2*U[j][i] + U[j][i-1]));
        }
    }    

    // Update U, FU, u, p, and c
    updateflow(N, Unew, U, FU, u, c, p);

}
/* ------------------------------ Lax-Wendroff ------------------------------ */
// Inputs: N = number of grid points, U = vector of conserved variables where
// U[0][:] = rho, U[1][:] = rho*u, U[2][:] = rho*E, FU = vector of flux vectors
// where FU[0][:] = rho*u, F[1][:] = rho*u^2+p, FU[2][:] = rho*u*H,
// u = primitive velocity, c = speed of sound, dt = delta time, dx = delta x.
// Outputs: Updates U, FU, u, p, and c for current timestep.

void laxwendroff(int N, double U[3][N], double FU[3][N], double u[N], double p[N],
                double c[N], double dt, double dx){ 
    // Initialize temporary U_new matrix, residual, aand flux
    double Unew[3][N];
    double residual[3][N];
    double flux[3][N-1];

    // Initialize residual and Unew
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < 3; j ++){
            residual[j][i] = 0.0; // Residual actually N columns
            flux[j][i] = 0.0;
        }
    }

    // Evaluate flux
    flux_lw(N, flux, U, FU, u, c, residual);

    // Loop over points to update
    for (int i = 1; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            Unew[j][i] = U[j][i] - dt/dx * residual[j][i];
        }
    }    

    // Enforce boundary conditions
    // NOTE: No boundary conditions, endpoints are constant in time

    // Update U, u, p, and c
    updateflow(N, Unew, U, FU, u, c, p);

}

/* -------------------------------- Van Leer -------------------------------- */
void vanleer(int N, double U[3][N], double FU[3][N], double u[N], double p[N],
            double c[N], double dt, double dx){
    // Initialize temporary U_new matrix, residual, aand flux
    double Unew[3][N];
    double residual[3][N];
    double flux[3][N-1];

    // Initialize residual and Unew
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < 3; j ++){
            residual[j][i] = 0.0; // Residual actually N columns
            flux[j][i] = 0.0;
        }
    }

    // Evaluate flux
    flux_vl(N, flux, U, FU, u, c, residual);

    // Loop over points to update
    for (int i = 1; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            Unew[j][i] = U[j][i] - dt/dx * residual[j][i];
        }
    }    

    // Enforce boundary conditions
    // NOTE: No boundary conditions, endpoints are constant in time

    // Update U, u, p, and c
    updateflow(N, Unew, U, FU, u, c, p);

} 

/* ------------------------------- MacCormack ------------------------------- */
void maccormack(int N, double U[3][N], double FU[3][N], double u[N], double p[N], 
                double c[N], double dt, double dx){
    // Initialize temporary U_new matrix, residual, aand flux
    double Unew[3][N];

    // MacCormack predictor-corrector
    double Upre[3][N-1];
    double FUpre[3][N-1];
    double Ucor[3][N-1];

    // Predictor step, find U bar
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            Upre[j][i] = U[j][i] - dt/dx * (FU[j][i+1] - FU[j][i]);
        }
    }
    // Update predicted u, p, c to find f(U) bar
    for (int i = 1; i < N-1; i++){
        u[i] = Upre[1][i] / Upre[0][i];
        p[i] = perfgas_p(Upre[2][i]/Upre[0][i], Upre[0][i], u[i]); 
        c[i] = sqrt(GAMMA * p[i] / Upre[0][i]);
    }
    // Find predicted f(U) bar
    for (int i = 0; i < N-1; i++){
        FUpre[0][i] = Upre[1][i];
        FUpre[1][i] = Upre[1][i] * u[i] + p[i];
        FUpre[2][i] = Upre[1][i] * (Upre[2][i] + p[i])/Upre[0][i];
    }

    // Corrector step, find U double bar
    for (int i = 1; i < N; i++){
        for (int j = 0; j < 3; j++){
            Ucor[j][i] = U[j][i] - dt/dx * (FUpre[j][i] - FUpre[j][i-1]);
        }
    } 

    // Loop over points to update
    for (int i = 1; i < N-1; i++){
        for (int j = 0; j < 3; j++){
            Unew[j][i] = 0.5 * (Upre[j][i] + Ucor[j][i]);
        }
    }    

    // Enforce boundary conditions
    // NOTE: No boundary conditions, endpoints are constant in time

    // Update U, u, p, and c
    updateflow(N, Unew, U, FU, u, c, p);
}