// header guards
#ifndef EULER_SOLVER_H
#define EULER_SOLVER_H

#include <vector>
#include "Variables.h"

class EulerSolver {
  private:
    int nx;
    double dx;
    double gamma;
    std::vector<ConservedVector> U;
   
  private:
    PrimitiveVector conservedToPrimitive(const ConservedVector& UC) const;
    ConservedVector primitiveToConserved(const PrimitiveVector& UP) const;
    ConservedVector evaluateFluxVector(const ConservedVector& FU) const;
    

  public:
    EulerSolver(int num_cells, double length); // constructor
    void initialiseState();
    void runSimulation(double t_end);

};

#endif