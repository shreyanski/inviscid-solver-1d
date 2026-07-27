#ifndef VARIABLES_H
#define VARIABLES_H

// data structure to hold primitive variables
struct PrimitiveVector {
  double rho  = 1.225;
  double u    = 0.0;
  double P    = 1.01325; // pressure
};

// data structure to hold the conserved variables
struct ConservedVector {
    double rho  = 1.0;   
    double rhou = 0.0;  
    double E    = 2.5; // energy which is f(rho,u,P,gamma)

    ConservedVector operator+(ConservedVector const& other) const {
      return {rho + other.rho, rhou + other.rhou, E + other.E};
    }
    ConservedVector operator-(ConservedVector const& other) const {
      return {rho - other.rho, rhou - other.rhou, E - other.E};
    }
    ConservedVector operator*(double scalar) const {
      return {rho * scalar, rhou * scalar, E * scalar};
    }

};


#endif // VARIABLES_H