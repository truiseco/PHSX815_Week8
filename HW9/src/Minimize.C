// C++ Standard Library Includes
#include <iostream>

// ROOT Includes
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// Directives and declarations for namespaces and namespace members
using std::exp, std::sin;

// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
double GausSin(const double *x){
  const double t = x[0];
  return exp(-t*t)*sin(t);
}

int main(int argc, char** argv){

   int randomSeed = -1;
   ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer();

   // set tolerance , etc...
   minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
   minimum->SetMaxIterations(10000);  // for GSL
   minimum->SetTolerance(0.001);
   minimum->SetPrintLevel(1);

   ROOT::Math::Functor f(&GausSin,1);
   double step[1] = {0.01};

   double variable[1] = { 0};
   if (randomSeed >= 0) {
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-20,20);
   }

   minimum->SetFunction(f);
   minimum->SetVariable(0,"x",variable[0], step[0]);
   minimum->Minimize();

   const double *xs = minimum->X();
   std::cout << "Minimum: f(" << xs[0] << "): "
             << minimum->MinValue()  << std::endl;

   if ( minimum->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
      std::cout << "Minimizer converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer failed to converge !!!" << std::endl;
      Error("NumericalMinimization","fail to converge");
   }

   return 0;
}
