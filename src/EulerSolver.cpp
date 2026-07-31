#include <iostream>
#include <cmath>
#include <functional>
#include "../include/EulerSolver.h"
#include "../include/Variables.h"

// Allocate memory & instantiate EulerSolver class (using constructor)
EulerSolver::EulerSolver(int num_cells, double length) 
: nx(num_cells), dx(length / num_cells), gamma(1.4) {
  U.resize(nx);
  F.resize(nx+1); // each cell has 2 shared faces except the boundary cells which have 1 shared face each
  
} 

// Sod-shock problem initialisation of primitive variables
void EulerSolver::initialiseState() {
  // initialise the conserved variable vector (using the primitive)
  for (int i=0; i<nx-1; ++i) {
    if (i < nx/2) {
      PrimitiveVector prim_left{1.225, 0.0, 1.0};
      U[i] = primitiveToConserved(prim_left);
    } else {
      PrimitiveVector prim_right{0.125, 0.0, 0.1};
      U[i] = primitiveToConserved(prim_right);
    }

    // debug
    std::cout << "Cell "<< i << " density: " << "rho = " << U[i].rho << std::endl;
  
  }
  
};

// Simulation execution method
void EulerSolver::runSimulation(double t_end, double CFL) {
  double t, dt {0.0};

  while (t < t_end) {

    double a, a_max, u {0.0};
   
    for (int i=0; i<nx; i++) {
      a = computeWaveSpeed(U[i]);
      u = conservedToPrimitive(U[i]).u;

      if (a + std::fabs(u) > a_max) {
        a_max = a + std::fabs(u);
      }
    }
    dt = CFL * dx / a_max;

    
    // applying MUSCL scheme to interpolate conserved variables at cell centers to faces
    
    // use low-order scheme at boundary cells [0] and [nx-1]
    Ufaces[0][0] = U[0];
    Ufaces[0][1] = U[0];
    Ufaces[nx-1][0] = U[nx-1];
    Ufaces[nx-1][1] = U[nx-1];
    
    ConservedVector epsilon {1e-8};
    // use higher-order scheme for interior faces
    for (int i=1; i<nx-1; i++) {
      // smoothness indicator calculation
      auto du_i_plus_half = U[i+1] - U[i]; // return type should be ConservedVector (contains rho, rhou, E)
      auto du_i_minus_half = U[i] - U[i-1]; 
      auto r_left = du_i_minus_half / (du_i_plus_half + epsilon);
      auto r_right = du_i_plus_half / (du_i_minus_half + epsilon); // return type is ConservedVector still? 

      // apply van Leer limiter
      ConservedVector r_left_fabs = fabs(r_left);
      ConservedVector psi_left = (r_left + r_left_fabs) / (ConservedVector{1.0,1.0,1.0} + r_left_fabs);
      ConservedVector r_right_fabs = fabs(r_right);
      ConservedVector psi_right = (r_right + r_right_fabs) / (ConservedVector{1.0,1.0,1.0} + r_right_fabs);

      Ufaces[i][0] = U[i] - 0.5 * psi_left * du_i_minus_half;
      Ufaces[i][1] = U[i] + 0.5 * psi_right * du_i_plus_half;

      // each interior shared face currently still has two U-values -> Rusanov flux for left & right faces?
      ConservedVector flux_left = Ufaces[i-1][0];
      ConservedVector flux_right = Ufaces[i][1];

      auto rho_left = flux_left.rho;
      auto rho_right = flux_right.rho;
      auto u_left = flux_left.rhou / rho_left;
      auto u_right = flux_right.rhou / rho_right;
      auto P_left = conservedToPrimitive(flux_left).P;
      auto P_right = conservedToPrimitive(flux_right).P;
      auto E_left = flux_left.E;
      auto E_right = flux_right.E;

      
      flux_left = {
        rho_left * u_left,
        P_left + rho_left * u_left * u_left,
        u_left * (E_left + P_left)
      };
      flux_right = {
        rho_right * u_right,
        P_right + rho_right * u_right * u_right,
        u_right * (E_right + P_right)
      };

      

    }


    t+=dt;
  }
 

};

// debug size and type of U
void EulerSolver::debug_grid() {
    
  // std::cout << Ufaces << std::endl;
    
}







//////////////////////////////////////////////////////////////////////////////////////
// PHYSICS/AUXILIARY METHODS IMPLEMENTED HERE
//////////////////////////////////////////////////////////////////////////////////////


PrimitiveVector EulerSolver::conservedToPrimitive(const ConservedVector& cons) const {
  PrimitiveVector prim;
  
  prim.rho = cons.rho;
  prim.u = cons.rhou / cons.rho;
  double kinetic_energy = 0.5 * prim.rho * prim.u * prim.u;
  prim.P = (gamma - 1.0) * (cons.E - kinetic_energy);

  return prim;
}

ConservedVector EulerSolver::primitiveToConserved(const PrimitiveVector& prim) const {
  ConservedVector cons;

  cons.rho  = prim.rho;
  cons.rhou = prim.rho * prim.u;
  double kinetic_energy = 0.5 * prim.rho * prim.u * prim.u;
  cons.E = prim.P / (gamma - 1.0) + kinetic_energy;

  return cons;
}

ConservedVector EulerSolver::evaluateCellFlux(const ConservedVector& cons) const {
  ConservedVector flux;
  PrimitiveVector prim = conservedToPrimitive(cons);

  flux.rho = cons.rho;
  flux.rhou = (cons.rhou * prim.u) + prim.P;
  flux.E = prim.u * (cons.E + prim.P);
  return flux; 
}

double EulerSolver::computeWaveSpeed(const ConservedVector& cons)  const {
  PrimitiveVector prim = conservedToPrimitive(cons);
  // only need speed of sound value in flux interface solver method (just return speed value here)
  return std::sqrt(gamma * prim.P / cons.rho);
}

// ConservedVector EulerSolver::evaluateFaceFlux(const ConservedVector& cons_left, const ConservedVector& cons_right) const {
//   ConservedVector faceflux;

  
// }