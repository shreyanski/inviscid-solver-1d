#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include "../include/EulerSolver.h"
#include "../include/Variables.h"
#include "../include/FluxSolver.h"

// solver class constructor and initialiser
EulerSolver::EulerSolver(int num_cells, double length) 
: nx(num_cells), 
  dx(length / num_cells),
  gamma(1.4) ,
  U(nx),
  Ufaces(nx), // coupled face indices and cell indices (nx cells means nx+1 faces)
  Ffaces(nx)
{ 
};

// Sod-shock problem initialisation of primitive variables
void EulerSolver::initialiseState() {
  
  for (int i=0; i<nx; ++i) {
    if (i < nx/2) {
      PrimitiveVector prim_left{1.225, 0.0, 1.0};
      U[i] = primitiveToConserved(prim_left);
    } else {
      PrimitiveVector prim_right{0.125, 0.0, 0.1};
      U[i] = primitiveToConserved(prim_right);
    }

  }
  
};

// std::tuple<ConservedVector, ConservedVector> EulerSolver::FluxSolver::interpolate_MUSCL(EulerSolver& solv) {
//   int nx = solv.nx;
//   double dx = solv.dx;
//   double gamma = solv.gamma;
//   std::vector<ConservedVector> U = solv.U;
//   std::vector<std::array<ConservedVector,2>> Ufaces = solv.Ufaces;
//   std::vector<std::array<ConservedVector,2>> Ffaces = solv.Ffaces;
  

//     // APPLYING MUSCL INTERPOLATION SCHEME 
//     // use low-order scheme for faces at boundary cells [0] and [nx-1]
//     Ufaces[0][0] = U[0]; // second index references selected face (left or right face)
//     Ufaces[0][1] = U[0]; 
//     Ufaces[nx-1][0] = U[nx-1];
//     Ufaces[nx-1][1] = U[nx-1];
//     ConservedVector epsilon {1e-8, 1e-8, 1e-8};
//     // use higher-order scheme for interior grid faces
//     for (int i=1; i<nx-1; i++) { // loops over interior cells
//       // calculate smoothness indicator 
//       auto du_i_plus_half = U[i+1] - U[i]; // return type is ConservedVector (contains rho, rhou, E variables)
//       auto du_i_minus_half = U[i] - U[i-1]; 
//       auto rL = du_i_minus_half / (du_i_plus_half + epsilon);
//       auto rR = du_i_plus_half / (du_i_minus_half + epsilon); 

//       // apply van Leer limiter
//       ConservedVector psiL = (rL + fabs(rL)) / (ConservedVector{1.0,1.0,1.0} + fabs(rL));
//       ConservedVector psiR = (rR + fabs(rR)) / (ConservedVector{1.0,1.0,1.0} + fabs(rR));

//       // conserved variables at faces for each cell i
//       Ufaces[i][0] = U[i] - 0.5 * psiL * du_i_minus_half;
//       Ufaces[i][1] = U[i] + 0.5 * psiR * du_i_plus_half;

//       // // DEBUGGING
//       // std::cout << "Ufaces_WEST_of_PREV: " << "ITER_" + std::to_string(iter) + 
//       // "_xi_" + std::to_string(i) << " --> " + 
//       // std::to_string(Ufaces[i][0].rho) << "\t" << std::to_string(Ufaces[i][0].rhou)
//       // << "\t" << std::to_string(Ufaces[i][0].E) << std::endl;

//     } // end loop over interior cells
    
//     // fluxes at boundary faces (leftmost & rightmost) -> compute from conserved variables at faces
//     ConservedVector qL_bound = Ufaces[0][0];
//     auto primL_bound = solv.conservedToPrimitive(qL_bound);
//     Ffaces[0][0] = {qL_bound.rhou,
//                     primL_bound.P + qL_bound.rhou * primL_bound.u, 
//                     primL_bound.u * (qL_bound.E + primL_bound.P)};
//     ConservedVector qR_bound = Ufaces[nx-1][1];
//     auto primR_bound = solv.conservedToPrimitive(qR_bound);
//     Ffaces[nx-1][1] = {qR_bound.rhou,
//                        primR_bound.P + qR_bound.rhou * primR_bound.u,
//                        primR_bound.u * (qR_bound.E + primR_bound.P)};

//     // each interior face currently has two U-values associated (face is shared w/ two cells)
//     for (int i=1; i<=nx-1; i++) { // loops over interior 

//       ConservedVector qL = Ufaces[i-1][1]; // right face of previous cell
//       ConservedVector qR = Ufaces[i][0];   // left face of current cell

//       auto rhoL = qL.rho;
//       auto rhoR = qR.rho;
//       auto uL = qL.rhou / rhoL;
//       auto uR = qR.rhou / rhoR;
//       auto PL = solv.conservedToPrimitive(qL).P;
//       auto PR = solv.conservedToPrimitive(qR).P;
//       auto EL = qL.E;
//       auto ER = qR.E;
//       auto aL = solv.computeWaveSpeedLocal(qL);
//       auto aR = solv.computeWaveSpeedLocal(qR);
    

//       ConservedVector fluxL = {
//         rhoL * uL,
//         PL + rhoL * uL * uL,
//         uL * (EL + PL)
//       };
//       ConservedVector fluxR = {
//         rhoR * uR,
//         PR+ rhoR * uR * uR,
//         uR * (ER + PR)
//       };

//       return {fluxL, fluxR, qL, qR};
// };


// Simulation execution method
void EulerSolver::runSimulation(double t_end, double CFL) {
  double t{0.0}, dt{0.0};
  int iter{0};

  while (t < t_end) {
    
    // ALLOWABLE TIMESTEP BASED ON CFL AND JACOBIAN
    auto Uold = U;
    double u{0.0}, a{0.0}, a_max{0.0};
    for (int i=0; i<nx; i++) {
      a = computeWaveSpeedLocal(U[i]);
      u = conservedToPrimitive(U[i]).u;

      if (a + std::fabs(u) > a_max) {
        a_max = a + std::fabs(u);
      }
    }
    dt = CFL * dx / a_max;

    // FluxSolver flux;
    // [ConservedVector fluxL, ConservedVector fluxR, ConservedVector qL, ConservedVector qR] = flux.interpolate_MUSCL(*this);

    // APPLYING MUSCL INTERPOLATION SCHEME 
    // use low-order scheme for faces at boundary cells [0] and [nx-1]
    Ufaces[0][0] = U[0]; // second index references selected face (left or right face)
    Ufaces[0][1] = U[0]; 
    Ufaces[nx-1][0] = U[nx-1];
    Ufaces[nx-1][1] = U[nx-1];
    ConservedVector epsilon {1e-8, 1e-8, 1e-8};
    // use higher-order scheme for interior grid faces
    for (int i=1; i<nx-1; i++) { // loops over interior cells
      // calculate smoothness indicator 
      auto du_i_plus_half = U[i+1] - U[i]; // return type is ConservedVector (contains rho, rhou, E variables)
      auto du_i_minus_half = U[i] - U[i-1]; 
      auto rL = du_i_minus_half / (du_i_plus_half + epsilon);
      auto rR = du_i_plus_half / (du_i_minus_half + epsilon); 

      // apply van Leer limiter
      ConservedVector psiL = (rL + fabs(rL)) / (ConservedVector{1.0,1.0,1.0} + fabs(rL));
      ConservedVector psiR = (rR + fabs(rR)) / (ConservedVector{1.0,1.0,1.0} + fabs(rR));

      // conserved variables at faces for each cell i
      Ufaces[i][0] = U[i] - 0.5 * psiL * du_i_minus_half;
      Ufaces[i][1] = U[i] + 0.5 * psiR * du_i_plus_half;

      // DEBUGGING
      std::cout << "Ufaces_WEST_of_PREV: " << "ITER_" + std::to_string(iter) + 
      "_xi_" + std::to_string(i) << " --> " + 
      std::to_string(Ufaces[i][0].rho) << "\t" << std::to_string(Ufaces[i][0].rhou)
      << "\t" << std::to_string(Ufaces[i][0].E) << std::endl;

    } // end loop over interior cells
    
    // fluxes at boundary faces (leftmost & rightmost) -> compute from conserved variables at faces
    ConservedVector qL_bound = Ufaces[0][0];
    auto primL_bound = conservedToPrimitive(qL_bound);
    Ffaces[0][0] = {qL_bound.rhou,
                    primL_bound.P + qL_bound.rhou * primL_bound.u, 
                    primL_bound.u * (qL_bound.E + primL_bound.P)};
    ConservedVector qR_bound = Ufaces[nx-1][1];
    auto primR_bound = conservedToPrimitive(qR_bound);
    Ffaces[nx-1][1] = {qR_bound.rhou,
                       primR_bound.P + qR_bound.rhou * primR_bound.u,
                       primR_bound.u * (qR_bound.E + primR_bound.P)};

    // each interior face currently has two U-values associated (face is shared w/ two cells)
    for (int i=1; i<=nx-1; i++) { // loops over interior 

      ConservedVector qL = Ufaces[i-1][1]; // right face of previous cell
      ConservedVector qR = Ufaces[i][0];   // left face of current cell

      auto rhoL = qL.rho;
      auto rhoR = qR.rho;
      auto uL = qL.rhou / rhoL;
      auto uR = qR.rhou / rhoR;
      auto PL = conservedToPrimitive(qL).P;
      auto PR = conservedToPrimitive(qR).P;
      auto EL = qL.E;
      auto ER = qR.E;
      auto aL = computeWaveSpeedLocal(qL);
      auto aR = computeWaveSpeedLocal(qR);
    

      ConservedVector fluxL = {
        rhoL * uL,
        PL + rhoL * uL * uL,
        uL * (EL + PL)
      };
      ConservedVector fluxR = {
        rhoR * uR,
        PR+ rhoR * uR * uR,
        uR * (ER + PR)
      };

      // DEBUGGING
      std::cout << "Ufaces_WEST_PREV: " << "ITER_" + std::to_string(iter) + 
      "_xi_" + std::to_string(i) << " --> ( " +
      std::to_string(Ufaces[i][1].rho) + " , " + std::to_string(Ufaces[i][1].rhou)
      + " , " + std::to_string(Ufaces[0][1].E) + " )" << std::endl;

      // need to consolidate these two fluxes into a single flux (for each face i)
      auto S_max = std::max(std::fabs(uL) + aL, std::fabs(uR) + aR);

      // // DEBUGGING
      // std::cout << "(global) current S_max: " << "ITER_" + std::to_string(iter) + 
      // "_xi_" + std::to_string(i) +
      // " --> " + std::to_string(S_max) << std::endl;

      ConservedVector flux_rusanov = 0.5 * (fluxL + fluxR) - 0.5 * S_max * (qR - qL);
      Ffaces[i-1][1] = flux_rusanov; // right face of cell i-1
      Ffaces[i][0] = flux_rusanov; // left face of cell i

      // DEBUGGING
      std::cout << "Ffaces_DENS_EAST_of_prev: ITER_" + std::to_string(iter) + 
      "_xi_" + std::to_string(i) + " --> " +
      std::to_string(Ffaces[i-1][1].rho) << std::endl;
      
    } // end for loop 


    
            
  for (int i=0; i<nx; i++) {
    U[i] = Uold[i] - (dt/dx) * (Ffaces[i][1] - Ffaces[i][0]);
  }
  // update boundary conditions
  U[0] = U[1];
  U[nx-1] = U[nx-2];
  
  t+=dt;
  iter++; 


  // write to file for csv post-processing
  writePrimitiveVariables(U, iter, t, t_end);



  //  // Print active tracking diagnostics cleanly to the console terminal
  //   std::cout << "Time: " << std::scientific << std::setprecision(3) << t 
  //             << " | Step: " << std::fixed << std::setprecision(5) << dt 
  //             << " | Max Wave Speed: " << a_max << "\r" << std::flush;
  

  } // end while loop

   
};





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

double EulerSolver::computeWaveSpeedLocal(const ConservedVector& cons)  const {
  PrimitiveVector prim = conservedToPrimitive(cons);
  // only need speed of sound value in flux interface solver method (just return speed value here)
  return std::sqrt(gamma * prim.P / cons.rho);
}

// ConservedVector EulerSolver::evaluateFaceFlux(const ConservedVector& cons_left, const ConservedVector& cons_right) const {
//   ConservedVector faceflux;

  
// }

// finally write the solution to .csv files
void EulerSolver::writePrimitiveVariables(const std::vector<ConservedVector>& prim, int iter, double t, double t_end) {
    // write to file for csv post-processing
    std::ofstream outputFile;
    // Only save structural CSV profiles occasionally (e.g., every 50 steps)
    if (iter % 50 == 0 || t >= t_end) {
      std::ofstream outputFile("build/Primitive_Iter_" + std::to_string(iter) + ".csv");
      outputFile << "x,rho,u,P\n"; // Using short strings over slow std::endl
      for (int i = 0; i < nx; i++) {
          auto prim = conservedToPrimitive(U[i]);
          outputFile << i * dx << "," << U[i].rho << "," << prim.u << "," << prim.P << "\n";
      }
      outputFile.close();
    }
  }

