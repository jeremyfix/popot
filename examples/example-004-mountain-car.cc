#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::algorithm::ParticleStochasticSPSO::VECTOR_TYPE TVector;

typedef popot::problems::rl::mountain_car::Params              StateParams;
typedef popot::problems::rl::mountain_car::State<StateParams>  State;
typedef popot::problems::rl::mountain_car::Action              Action;


typedef popot::problems::rl::mountain_car::RBFParams RBFParams;
typedef popot::problems::rl::mountain_car::RBFPolicy<StateParams, RBFParams> Policy;

typedef popot::problems::rl::mountain_car::Simulator<Policy, StateParams>  Simulator;

#define NB_MONTECARLO 30
#define MAX_EPOCH 50
#define NB_POLICY_EVALUATIONS 1000

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  size_t dimension = Policy::dimension();

  auto stop =   [] (double fitness, int epoch) -> bool { return epoch >= MAX_EPOCH;};
  auto cost_function = [dimension] (TVector& pos) -> double { 
    return -Simulator::evaluate(pos.getValuesPtr());
  };

  auto algo = popot::algorithm::stochastic_montecarlo_spso2006(Policy::dimension(), 
							       Policy::lbound, 
							       Policy::ubound, 
							       stop, 
							       cost_function, 
							       NB_MONTECARLO);
    
  // We run the algorithm
  algo->run(1);

  std::cout << "Best particle :" << algo->getBest() << std::endl;

  std::cout << "Mean trajectory length over " << NB_POLICY_EVALUATIONS << " runs : " 
	    << Simulator::test_policy(algo->getBest().getPosition().getValuesPtr(), NB_POLICY_EVALUATIONS) << std::endl;
  std::cout << "The episode ending at " << StateParams::max_length_episode() << " steps " << std::endl;
    

  delete algo;
}
