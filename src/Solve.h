#ifndef Solve_H_
#define Solve_H_

const double GAMMA; // Specific heat capacity
const double RAIR;  // Specific gas constant, air

/* ------------------------------- Perfect Gas ------------------------------ */
double perfgas_p(double E, 
                double rho,
                double u);

double perfgas_e(double p,
                double rho,
                double u);

/* ------------------------------ CFL Condition ----------------------------- */
double dtcalc(double cfl, double dx, double umax);

/* --------------------------- Boundary Conditions -------------------------- */
double bc_none(int N, double U[3][N], double FU[3][N]);

/* ----------------------------- Flux Functions ----------------------------- */
void flux_lw(int N,
            double flux[3][N-1], 
            double U[3][N], 
            double FU[3][N], 
            double u[N], 
            double c[N],
            double residual[3][N]);

void flux_vl(int N,
            double flux[3][N-1], 
            double U[3][N], 
            double FU[3][N], 
            double u[N], 
            double c[N],
            double residual[3][N]);

/* ---------------------------- Update Flow Field --------------------------- */
void updateflow(int N,
                double Unew[3][N],
                double U[3][N],
                double FU[3][N], 
                double u[N],
                double c[N],
                double p[N]);

/* -------------------------------------------------------------------------- */
/* ----------------------- Compressible Euler Solvers ----------------------- */
/* -------------------------------------------------------------------------- */

/* --------------------------- Central Difference --------------------------- */
void centraldiff(int N,
                double U[3][N],
                double FU[3][N],
                double u[N],
                double p[N],
                double c[N],
                double dt,
                double dx);

/* ------------------------------ Lax-Wendroff ------------------------------ */
void laxwendroff(int N,
                double U[3][N],
                double FU[3][N],
                double u[N],
                double p[N],
                double c[N],
                double dt,
                double dx);

/* -------------------------------- Van Leer -------------------------------- */
void vanleer(int N,
            double U[3][N],
            double FU[3][N],
            double u[N],
            double p[N],
            double c[N],
            double dt,
            double dx);

/* ------------------------------- MacCormack ------------------------------- */
void maccormack(int N,
            double U[3][N],
            double FU[3][N],
            double u[N],
            double p[N],
            double c[N],
            double dt,
            double dx);

#endif