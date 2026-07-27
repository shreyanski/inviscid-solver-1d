#include <iostream>
#include "../include/Variables.h"
#include "../include/EulerSolver.h"

int main() {
  
  // ConservedVector cell_1{2.5, 5.0, 10.0};
  // ConservedVector cell_2{1.5, 3.0, 6.0};
  // ConservedVector cell_sum = cell_1 + cell_2;
  // ConservedVector cell_sub = cell_1 - cell_2;
  // ConservedVector cell_scal = cell_1 * 10.0;
  // std::cout << "Testing operator overloading for ConservedVector struct" << std::endl;
  // std::cout << "\t Sum of cell densities: " << cell_sum.rho <<  " kg/m^3" << std::endl;
  // std::cout << "\t Diff of cell energy densities: " << cell_sub.E << " J/m^3" << std::endl;
  // std::cout << "\t Scaled cell momentum density: " << cell_scal.rhou << " kg/(m^2*s)" << std::endl;

  std::cout << "========= CFD SOLVER SANITY CHECK============ " << std::endl;
  // Construct the solver
  int test_cells = 5;
  double test_length = 1.0;
  EulerSolver test_solver(test_cells, test_length);

  // Run the initialise member function
  std::cout << "Initialising the simulation grid... " << std::endl;
  test_solver.initialiseState();

  // Run the simulation 
  test_solver.runSimulation(1.0);
  std::cout << "==========  CFD WORKING IF YOU GOT HERE!  ===========" << std::endl;

  
  return 0;
}