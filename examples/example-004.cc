#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::problems::Ackley Problem;

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  size_t dimension = 30;
  Problem p(dimension);

  auto lbound = [&p] (size_t index) -> double { return p.get_lbound(index); };
  auto ubound = [&p] (size_t index) -> double { return p.get_ubound(index); };
  auto stop =   [&p] (double fitness, int epoch) -> bool { return p.stop(fitness, epoch);};
  auto cost_function = [&p] (double* pos) -> double { return p.evaluate(pos);};

  auto algo = popot::algorithm::harmony(dimension, lbound, ubound, stop, cost_function);

  algo->run();

  delete algo;
}
