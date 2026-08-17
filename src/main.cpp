#include <iostream>
#include "../include/Variables.h"
#include "../include/EulerSolver.h"

int main() {


  // Construct the solver
  std::cout << "========= CFD SOLVER SANITY CHECK============ " << std::endl;
  int test_cells = 100;
  double test_length = 1.0;
  double test_time = 0.2;
  double test_CFL = 0.1;

  std::cout << "--> Allocating memory for the grid of conserved vectors..." << std::endl;
  EulerSolver debugger(test_cells, test_length);
  std::cout << "--> Allocation successful!" << std::endl;


  std::cout << "--> Initializing the simulation grid... " << std::endl;
  debugger.initialiseState();
  std::cout << "--> Initialization successful!" << std::endl;


  // debugger.debug_U(test_cells);

  std::cout << "--> Running the simulation..." << std::endl;
  debugger.runSimulation(test_time, test_CFL);

  std::cout << "==========  SUCCESS!!!  ===========\n";
  return 0;
}