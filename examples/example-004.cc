#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::algorithm::ParticleStochasticSPSO::VECTOR_TYPE TVector;
//typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
typedef popot::problems::QuarticNoise Problem;
//typedef popot::problems::Rosenbrock Problem;

int main(int argc, char* argv[]) {



  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  size_t dimension = 2;
  Problem p(dimension);

  auto lbound = [&p] (size_t index) -> double { return p.get_lbound(index); };
  auto ubound = [&p] (size_t index) -> double { return p.get_ubound(index); };
  auto stop =   [&p] (double fitness, int epoch) -> bool { return p.stop(fitness, epoch);};
  auto cost_function = [&p] (TVector& pos) -> double { return p.evaluate(pos.getValuesPtr());};

  //auto algo = popot::algorithm::spso2006(15,dimension, lbound, ubound, stop, cost_function);
  auto algo = popot::algorithm::stochastic_spso2006(30,dimension, lbound, ubound, stop, cost_function);

  algo->run(1);

  std::cout << "Best particle :" << algo->getBest() << std::endl;

  delete algo;

  /*
  popot::PSO::particle::StochasticParticle<> part(dimension);
  popot::initializer::position::uniform_random(part.getPosition(), lbound, ubound);
  part.evaluateFitness(cost_function);
  part.evaluateFitness(cost_function);
  std::cout << part << std::endl;

  popot::PSO::particle::compareFitnessMonteCarlo(part, part, 10, cost_function);

  std::ofstream outfile("toto.save");
  part.save(outfile);
  outfile.close();

  
  std::ifstream infile("toto.save");
  popot::PSO::particle::StochasticParticle<> part2(dimension);
  part2.load(infile);
  infile.close();

  std::cout << part2 << std::endl;
  */
}
