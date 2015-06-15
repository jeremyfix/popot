#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::Vector<double> TVector;

int main(int argc, char * argv[]) {
  
  size_t dimension = 10;
  auto lbound = [](size_t dim) -> double { return -5;};
  auto ubound = [](size_t dim) -> double { return 5;};
  auto stop = [](double fitness, int epoch) -> bool { return fitness <= 1e-10;};
  auto cost = [dimension](TVector& pos) -> double {
    double res = 0.0;
    for(unsigned int i = 0 ; i < dimension ; ++i) {
      double x = pos.getValuesPtr()[i];
      res += (x-1)*(x-1);
    }
    return res;
  };

  size_t nbloups = 1000;
  size_t max_iter = 100;
  auto algo = popot::algorithm::gwo(nbloups, dimension, lbound, ubound,
				    max_iter, stop, cost);

  std::ofstream outfile("fitnesses.data");
  for(unsigned int i = 0 ; i < 100 ; ++i) {
    algo->step();
    algo->print(1);
    outfile << algo->getBestFitness() << std::endl;
  }
  outfile.close();
  delete algo;
}
