// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"



popot::problems::Base* make_problem(std::string problem_name, int dimension) {
  if(problem_name == "Ackley")
    return new popot::problems::Ackley(dimension);
  else if(problem_name == "Quadric")
    return new popot::problems::Quadric(dimension);
  else
    return new popot::problems::Griewank(dimension);
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();

  if(argc != 3) {
    std::cerr << "Usage : " << argv[0] << " problemname dimension " << std::endl;
    return -1;
  }

  std::string problem_name(argv[1]);
  int dimension = atoi(argv[2]);

  auto problem = make_problem(problem_name, dimension);
  std::cout << problem->get_ubound(0) << std::endl;
  
}
