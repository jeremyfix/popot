// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::Vector<double> TVector;

popot::problems::Base* make_problem(std::string problem_name, int dimension) {
  if(problem_name == "Ackley")
    return new popot::problems::Ackley(dimension);
  else if(problem_name == "Quadric")
    return new popot::problems::Quadric(dimension);
  else if(problem_name == "Griewank")
    return new popot::problems::Griewank(dimension);
  else if(problem_name == "Sphere")
    return new popot::problems::Sphere(dimension);
  else if(problem_name == "Rastrigin")
    return new popot::problems::Rastrigin(dimension);
  else if(problem_name == "Rosenbrock")
    return new popot::problems::Rosenbrock(dimension);
  else if(problem_name == "Schwefel1_2")
    return new popot::problems::Schwefel1_2(dimension);
  else if(problem_name == "Schwefel")
    return new popot::problems::Schwefel(dimension);
  else if(problem_name == "Salomon")
    return new popot::problems::Salomon(dimension);
  else if(problem_name == "Dropwave") {
    if(dimension != 2)
      throw std::logic_error("Dropwave is only in dimension 2");
    return new popot::problems::Dropwave();
  }
  else 
    throw std::logic_error("Unrecognized problem name");
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();

  if(argc != 4) {
    std::cerr << "Usage : " << argv[0] << " problemname dimension algorithm" << std::endl;
    std::cerr << "With : "
	      << " - problemname in {Ackley, Quadric, Griewank, Sphere, Rastrigin, Rosenbrock, Schwefel1_2, Schwefel, Salomon}"
	      << " - algorithm in {spso2011, spso2006, abc, gwo} "
	      << std::endl;
    return -1;
  }

  std::string problem_name(argv[1]);
  int dimension = atoi(argv[2]);
  std::string algo_name(argv[3]);

  auto problem = make_problem(problem_name, dimension);
  auto lbound = [problem] (size_t index) -> double { return problem->get_lbound(index); };
  auto ubound = [problem] (size_t index) -> double { return problem->get_ubound(index); };
  auto stop =   [problem] (double fitness, int epoch) -> bool { return problem->stop(fitness, epoch);};
  auto cost_function = [problem] (TVector &pos) -> double { return problem->evaluate(pos.getValuesPtr());};

  popot::algorithm::Base* algo;
  if(algo_name == "spso2011")
    algo = popot::algorithm::spso2011(dimension, lbound, ubound, stop, cost_function);
  else if(algo_name == "spso2006")
    algo = popot::algorithm::spso2006(dimension, lbound, ubound, stop, cost_function);
  else if(algo_name == "abc") 
    algo = popot::algorithm::abc(50, dimension, lbound, ubound, stop, cost_function);
  else if ( algo_name == "gwo"){
    std::cout << "On lance l'algo gwo " << std::endl; //test
    size_t nb_max_iterations = problem->max_fe();
    algo = popot::algorithm::gwo(1000, dimension, lbound, ubound, nb_max_iterations, stop, cost_function);

  }
  else {
    std::cerr << "Invalid algorithm name " << std::endl;
    return -1;
  }
  if(!algo) 
    throw std::logic_error("The algorithm is not initialized");

  auto benchmark = popot::benchmark::make_benchmark(*algo, *problem, 2, 1);
  benchmark.run(1);
  std::cout << benchmark << std::endl;

  std::string result_filename_base = problem_name + "_" + std::to_string(dimension) + "_" + algo_name;
  unsigned int index = 0;
  std::string result_filename = result_filename_base + "_" + std::to_string(index++) + ".json";
  std::ifstream infile(result_filename.c_str());
  while(infile.is_open()) {
    infile.close();
    result_filename = result_filename_base + "_" + std::to_string(index++) + ".json";
    infile.open(result_filename.c_str());
  }
  infile.close();


  benchmark.dump_json(result_filename, problem_name, algo_name);
  std::cout << "Results saved in " << result_filename << std::endl;

  
}
