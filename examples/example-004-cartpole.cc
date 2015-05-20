#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::algorithm::ParticleStochasticSPSO::VECTOR_TYPE TVector;

typedef popot::problems::rl::inverted_pendulum::Params              StateParams;
typedef popot::problems::rl::inverted_pendulum::Pendulum            PendulumParams;
typedef popot::problems::rl::inverted_pendulum::State<StateParams>  State;
typedef popot::problems::rl::inverted_pendulum::Action              Action;


//typedef popot::problems::rl::inverted_pendulum::NormalPolicy<State> Policy;

typedef popot::problems::rl::inverted_pendulum::RBFParams RBFParams;
typedef popot::problems::rl::inverted_pendulum::RBFPolicy<State, RBFParams> Policy;

typedef popot::problems::rl::inverted_pendulum::Simulator<Policy, StateParams, PendulumParams>  Simulator;

#define NB_MONTECARLO 1
#define MAX_EPOCH 100
#define NB_POLICY_EVALUATIONS 500

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
