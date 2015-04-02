// Example of how to make use of the ABC algorithm on a custom evaluate function


// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

// Define the vector type and the problem
typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
typedef popot::problems::Ackley Problem;

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
 
  // Some initialization of static fields
  size_t dimension = 30;
  Problem p(dimension);
  
  // Let's create a swarm 
  // we might use spso2006, spso2007 or spso2011

  auto lbound = [&p] (size_t index) -> double { return p.get_lbound(index); };
  auto ubound = [&p] (size_t index) -> double { return p.get_ubound(index); };
  auto stop =   [&p] (double fitness, int epoch) -> bool { return p.stop(fitness, epoch);};
  auto cost_function = [&p] (TVector &pos) -> double { return p.evaluate(pos.getValuesPtr());};


  size_t colony_size = 50; 
  auto algo = popot::algorithm::abc(colony_size, dimension,
				    lbound, ubound,
				    stop, cost_function);

  algo->init();
  algo->run();

  std::cout << "Best minimum found :" << algo->getBest().getFValue() << " in " << algo->getEpoch() << " steps " << std::endl;
  std::cout << "Position of the optimum : " << algo->getBest() << std::endl;
  std::cout << std::endl;
  std::cout << p.getFE() << " Function Evaluations" << std::endl;

  delete algo;
  
}
