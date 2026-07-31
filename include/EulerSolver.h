// header guards
#ifndef EULER_SOLVER_H
#define EULER_SOLVER_H

#include <vector>
#include <array>
#include "Variables.h"

class EulerSolver {
  private:
    int nx;
    double dx;
    double gamma;
    
    std::vector<ConservedVector> U;
    std::vector<std::array<ConservedVector, 2>> Ufaces;
    std::vector<ConservedVector> F; // use as flux vector 
   
  private:
    PrimitiveVector conservedToPrimitive(const ConservedVector& cons) const;
    ConservedVector primitiveToConserved(const PrimitiveVector& prim) const;
    ConservedVector evaluateCellFlux(const ConservedVector& cellflux) const;
    double computeWaveSpeed(const ConservedVector& cons) const;
    ConservedVector evaluateFaceFlux(const ConservedVector& cons_left, const ConservedVector& cons_right) const;



  public:
    EulerSolver(int num_cells, double length); // constructor
    void initialiseState();
    void runSimulation(double t_end, double CFL);

    void debug_grid();

};

#endif