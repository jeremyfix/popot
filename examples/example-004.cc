#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

//typedef popot::algorithm::StochasticParticleSPSO::VECTOR_TYPE TVector;
typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
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
  auto cost_function = [&p] (TVector& pos) -> double { return p.evaluate(pos.getValuesPtr());};

  //auto algo = popot::algorithm::stochastic_spso2006(dimension, lbound, ubound, stop, cost_function);
  auto algo = popot::algorithm::spso2006(dimension, lbound, ubound, stop, cost_function);

  algo->run(1);

  delete algo;
}
