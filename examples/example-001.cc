// In this example, we show a basic usage of the SPSO algorithms, using the default ones
// provided by the library

// We first define the generator of random numbers
// THIS is absolutely necessary to define RNG_GENERATOR before
// including popot.h as some of the codes require it to be defined
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define the vector type and the problem
typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
typedef popot::problems::Rosenbrock Problem;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  size_t dimension = 5;
  Problem p(5);
  
  // Let's create a swarm 
  // we might use spso2006, spso2007 or spso2011
  auto algo = popot::algorithm::spso2011(dimension,
  					 [&p] (size_t index) -> double { return p.get_lbound(index); },
  					 [&p] (size_t index) -> double { return p.get_ubound(index); },
  					 [&p] (double fitness, int epoch) -> bool { return p.stop(fitness, epoch);},
  					 [&p] (TVector &pos) -> double { return p.evaluate(pos.getValuesPtr());}
					 );

  // Let us save the particles
  algo.save("particles.data");
  
  // We now run our algorithm
  algo.run(1);

  // We load the particles we had initially
  algo.load("particles.data");

  std::cout << "Nb steps : " << algo.epoch << " ; # of generated neighborhoods : " << algo.nb_new_neigh << std::endl;

}
