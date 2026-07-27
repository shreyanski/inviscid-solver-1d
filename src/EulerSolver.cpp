#include <iostream>
#include "../include/EulerSolver.h"

EulerSolver::EulerSolver(int num_cells, double length) 
: nx(num_cells), dx(length / num_cells), gamma(1.4) {
  U.resize(nx);
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
    std::cout << "Cell "<< i <<" density values:" << std::endl;
    std::cout << "rho = " << U[i].rho << std::endl;
  }
  // debug
  std::cout << "Simulation grid initialised for the Sod-shock problem" << std::endl;

};

void EulerSolver::runSimulation(double t_end) {
  double t = 0.0;
  // while (t < t_end) {
  // debug
    std::cout << "Running time integration loop..." << std::endl;
  //   t+=1.0;
  // }
};


// PHYSICS/AUXILIARY METHODS IMPLEMENTED
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


