#ifndef VARIABLES_H
#define VARIABLES_H

#include <cmath>

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
   
    ConservedVector operator*(ConservedVector const& other) const {
      return {rho * other.rho, rhou * other.rhou, E * other.E};
    }
    ConservedVector operator/(ConservedVector const& other) const {
      return {rho / other.rho, rhou / other.rhou, E / other.E};
    } 


};

// allows inline operator overload for commutative multiplication btw ConservedVector and scalar (must be free functions, not member)
inline ConservedVector operator*(const ConservedVector& cons, double scalar) {
  return ConservedVector{
    cons.rho * scalar,
    cons.rhou * scalar,
    cons.E * scalar
  };
}
inline ConservedVector operator*(double scalar, const ConservedVector& cons) {
  return cons * scalar; // reuse the function above since its commutative
}

// allow fabs to accept user-defined data type (defining it in the same namespace as ConservedVector)
inline ConservedVector fabs(const ConservedVector& cons) {
  return ConservedVector{
    std::fabs(cons.rho),
    std::fabs(cons.rhou),
    std::fabs(cons.E)
  };
}


#endif // VARIABLES_H