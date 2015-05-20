#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::algorithm::ParticleStochasticSPSO::VECTOR_TYPE TVector;



int main(int argc, char* argv[]) {



  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  size_t dimension = 2;

  auto lbound = [] (size_t index) -> double { return -5;};
  auto ubound = [] (size_t index) -> double { return  5;};
  auto stop =   [] (double fitness, int epoch) -> bool { return epoch >= 1000;};
  auto cost_function = [dimension] (TVector& pos) -> double { 
    double fitness = 0.0;
    for(size_t i = 0 ; i < dimension ; ++i)
      fitness += (pos[i]-1.5) * (pos[i]-1.5);
    fitness += popot::math::normal(0.0, 1.0);
    return fitness;
  };

  auto algo = popot::algorithm::stochastic_montecarlo_spso2011(dimension, lbound, ubound, stop, cost_function, 100);

  // We save the initial positions
  std::ofstream outfile;
  outfile.open("init.data");
  for(auto& p: algo->getParticles()) {
    for(size_t i = 0 ; i < dimension ; ++i)
      outfile << p.getPosition()[i] << " ";
    outfile << std::endl;
  }
  outfile.close();
  std::cout << "Positions saved in init.data" << std::endl;
    
  // We run the algorithm
  algo->run(1);
  std::cout << "Best particle :" << algo->getBest() << std::endl;

  // And save the final positions
  outfile.open("final.data");
  for(auto& p: algo->getParticles()) {
    for(size_t i = 0 ; i < dimension ; ++i)
      outfile << p.getPosition()[i] << " ";
    outfile << std::endl;
  }
  outfile.close();
  std::cout << "Positions saved in final.data" << std::endl;

  std::cout << "You can plot the positions using gnuplot : " << std::endl;
  std::cout << "   plot 'init.data' using 1:2, 'final.data' using 1:2" << std::endl;

  delete algo;
}
