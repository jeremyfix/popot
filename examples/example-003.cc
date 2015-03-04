// Example of how to make use of the ABC algorithm on a custom evaluate function


// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

double evaluate(double * x, size_t dimension, size_t & count)
{
  double * params = (double*) x;
  double fit = 0.0;
  double y_i, y_i_1;
  for(size_t i = 0 ; i < dimension-1 ; ++i)
    {
      y_i = params[i];
      y_i_1 = params[i+1];
      fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
    }
  count ++;

  return fit;
}

bool stop(double f, int epoch)
{
  return (f <= 1e-2) || (epoch >= 10000);
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  size_t colony_size = 50;
  size_t dimension = 15;
  size_t count=0;
  
  auto algo = popot::algorithm::abc(colony_size, dimension,
				    [] (size_t index) -> double { return -10; },
				    [] (size_t index) -> double { return  10; },
				    stop,
				    [&count,dimension] (double * params) -> double { return evaluate(params, dimension, count);});

  algo.init();
  algo.run();

  std::cout << "Best minimum found :" << algo.getBest().getFValue() << " in " << algo.getEpoch() << " steps " << std::endl;
  std::cout << "Position of the optimum : " << algo.getBest() << std::endl;
  std::cout << std::endl;
  std::cout << (int)count << " Function Evaluations" << std::endl;
  
}
