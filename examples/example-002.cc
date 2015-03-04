// In this example, we show how ones can define its own problem with a PSO
// The problem we try to solve is optimizing a Multilayer Perceptron to detect if a number is odd or even
// The input is represented by 7 segments indexed the following way
//         ----- 0 -----
//         |           |
//         5           1
//         |           |
//         ----- 6 -----
//         |           |
//         4           2
//         |           |
//         ----- 3 -----
// Therefore, we code :
//   = [ 6 5 4 3 2 1 0]  <- Index of the segments
// 0 = [ 0 1 1 1 1 1 1] ; 0x3F
// 1 = [ 0 0 0 0 1 1 0] ; 0x06
// 2 = [ 1 0 1 1 0 1 1] ; 0x5B
// 3 = [ 1 0 0 1 1 1 1] ; 0x4F
// 4 = [ 1 1 0 0 1 1 0] ; 0x66
// 5 = [ 1 1 0 1 1 0 1] ; 0x6D
// 6 = [ 1 1 1 1 1 0 1] ; 0x7D
// 7 = [ 0 0 0 0 1 1 1] ; 0x07
// 8 = [ 1 1 1 1 1 1 1] ; 0x7F
// 9 = [ 1 1 0 1 1 1 1] ; 0x6F
// In the code below, these are stored in inputs_mlp in the order : 0 1 2 3 4 5 6 ; 0 1 2 3 4 5 6 ....

#include <stdio.h>
#include <cmath>
#include <list>

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"


// Define our problem
// We define it as a MLP tested on our dataset
class MLPClassifier
{
public:
  const size_t nb_inputs;
  const size_t nb_outputs;
  const size_t nb_hidden;
  const size_t nb_digits;

  const size_t nb_parameters;// = (nb_inputs+1)*nb_hidden + (nb_hidden+1)*nb_outputs;
  int *hexa_codes;
  int *inputs_mlp;

  MLPClassifier(void) : 
    nb_inputs(7), nb_outputs(1), nb_hidden(4), nb_digits(10),
    nb_parameters((nb_inputs+1)*nb_hidden + (nb_hidden+1)*nb_outputs)
  {
    hexa_codes = new int[nb_digits];
    hexa_codes[0] = 0x3F;
    hexa_codes[1] = 0x06;
    hexa_codes[2] = 0x5B;
    hexa_codes[3] = 0x4F;
    hexa_codes[4] = 0x66;
    hexa_codes[5] = 0x6D;
    hexa_codes[6] = 0x7D;
    hexa_codes[7] = 0x07;
    hexa_codes[8] = 0x7F;
    hexa_codes[9] = 0x6F;

    inputs_mlp = new int[nb_inputs*nb_digits];
    size_t y;
    for(size_t x = 0 ; x < nb_digits ; ++x)
      {
	y = hexa_codes[x];
	for(size_t i = 0 ; i < nb_inputs ; ++i)
	  {
	    inputs_mlp[nb_inputs*x + i] = (y & 0x01);
	    y = y >> 1;
	  }
      }
  }

  ~MLPClassifier(void)
  {
    delete[] hexa_codes;
    delete[] inputs_mlp;
  }

  double get_lbound(size_t index)
  {  
    return -10;
  }

  double get_ubound(size_t index)
  {    
    return 10;
  }

  void display_digits(std::vector<int> &list_digits)
  {
    if(list_digits.size() == 0)
      {
	printf("Nothing to display .. \n");
	return;
      }
    // We want to display several digits on a row under their 7 segments format
    // We therefore display one after the other the rows for all the digits

    // Horizontal segments span 13 characters
    // Vertical segments span 3 characters
    size_t span_horiz = 13;
    size_t span_vert = 3;
    size_t space_digit = 3;

    // Segment 0
    for(size_t i = 0 ; i < list_digits.size() ; ++i)
      {
	for(size_t j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+0])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(size_t j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }
    printf("\n");
    // Segments 5 and 1
    for(size_t i = 0 ; i < span_vert ; ++i)
      {
	for(size_t j = 0 ; j < list_digits.size() ; ++j)
	  {
	    if(inputs_mlp[list_digits[j]*nb_inputs+5])
	      printf("|");
	    else
	      printf(" ");
	    for(size_t k = 0 ; k < span_horiz - 2 ; ++k)
	      printf(" ");
	    if(inputs_mlp[list_digits[j]*nb_inputs+1])
	      printf("|");
	    else
	      printf(" ");
	    for(size_t k = 0 ; k < space_digit ; ++k)
	      printf(" ");
	  }
	printf("\n");
      }

    // Segment 6
    for(size_t i = 0 ; i < list_digits.size() ; ++i)
      {
	for(size_t j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+6])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(size_t j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }
    printf("\n");
    // Segments 4 2
    for(size_t i = 0 ; i < span_vert ; ++i)
      {
	for(size_t j = 0 ; j < list_digits.size() ; ++j)
	  {
	    if(inputs_mlp[list_digits[j]*nb_inputs+4])
	      printf("|");
	    else
	      printf(" ");
	    for(size_t k = 0 ; k < span_horiz - 2 ; ++k)
	      printf(" ");
	    if(inputs_mlp[list_digits[j]*nb_inputs+2])
	      printf("|");
	    else
	      printf(" ");
	    for(size_t k = 0 ; k < space_digit ; ++k)
	      printf(" ");
	  }
	printf("\n");
      }
    // Segment 3
    for(size_t i = 0 ; i < list_digits.size() ; ++i)
      {
	for(size_t j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+3])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(size_t j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }

    printf("\n");
  }

  double transfer_function(double x)
  {
    // We use a sigmoidal transfer function
    return 1.0 / (1.0 + exp(-x));
  }

  double compute_output(int input_index, double * params)
  {
    double act_input[nb_inputs];
    double act_hidden[nb_hidden];
    double act_output;

    // Set up the input
    for(size_t i = 0 ; i < nb_inputs ; ++i)
      act_input[i] = inputs_mlp[nb_inputs * input_index + i];

    // Compute the activity of the hidden nodes
    size_t params_index = 0;
    for(size_t i = 0 ; i < nb_hidden ; ++i)
      {
	act_hidden[i] = params[params_index++];
	for(size_t j = 0 ; j < nb_inputs ; ++j)
	  {
	    act_hidden[i] += params[params_index++] * act_input[j];
	  }
	// Apply the transfer function
	act_hidden[i] = transfer_function(act_hidden[i]);
      }
    // Compute the activity of the output node
    act_output = params[params_index++];
    for(size_t j = 0 ; j < nb_hidden ; ++j)
      act_output += params[params_index++] * act_hidden[j];
    act_output = transfer_function(act_output);

    return act_output;
  }

  bool stop(double fitness, int epoch)
  {
    return (fitness <= 1e-4) || (epoch >= 1000);
  }

  double evaluate(double * params)
  {
    double fitness = 0.0;
    double act_output;
    // Test over all the inputs
    for(size_t k = 0 ; k < nb_digits ; ++k)
      {
	act_output = compute_output(k, params);

	fitness += pow(act_output - (k%2),2.0);
      }
    return fitness;
  }
};

// Let's typedef the problem so that the following is identical to example-001
typedef MLPClassifier Problem;
typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();

  // Create an instance of the problem
  Problem p;

  // For testing
  std::vector<int> odd_numbers;
  std::vector<int> even_numbers;
  std::vector<int> unclassified_numbers;
  double * params = new double[p.nb_parameters];

  // Let's create our swarm
  std::cout << " Nb params :" << p.nb_parameters << std::endl;
  auto algo = popot::algorithm::spso2011(p.nb_parameters,
  					 [&p] (size_t index) -> double { return p.get_lbound(index); },
  					 [&p] (size_t index) -> double { return p.get_ubound(index); },
  					 [&p] (double fitness, int epoch) -> bool { return p.stop(fitness, epoch);},
  					 [&p] (TVector &pos) -> double { return p.evaluate(pos.getValuesPtr());}
					 );

  // Let's generate the graph of the connections within the swarm
  algo.generateGraph("connections.dot");


  ////////////////////////////////////////:
  // Test before learning :
  printf("--------------------------------------------- \n Before learning : \n");
  printf("Best fitness : %f \n", algo.getBest().getFitness());
  
  for(size_t i = 0 ; i < p.nb_parameters ; ++i)
    params[i] = algo.getBest().getPosition()[i];

  for(size_t i = 0 ; i < p.nb_digits; ++i)
    {
      if(p.compute_output(i, params) >= 0.85)
	odd_numbers.push_back(i);
      else if(p.compute_output(i, params) <= 0.15)
	even_numbers.push_back(i);
      else
	unclassified_numbers.push_back(i);
    }

  printf("I classified %i numbers as EVEN : \n", even_numbers.size());
  p.display_digits(even_numbers);
  printf("I classified %i numbers as ODD : \n", odd_numbers.size());
  p.display_digits(odd_numbers);
  printf("I was not able to classify %i numbers : \n", unclassified_numbers.size());
  p.display_digits(unclassified_numbers);  

  ////////////////////////////////////////:
  // We now iterate the algorithm
  algo.run();
  std::cout << "epoch : " << algo.epoch << std::endl;
  std::cout << "\n" << std::endl;

  ////////////////////////////////////////:
  // Test after learning :
  printf("--------------------------------------------- \n After learning : \n");
  printf("Best fitness : %f \n", algo.getBest().getFitness());

  // Get the best parameters
  for(size_t i = 0 ; i < p.nb_parameters ; ++i)
    params[i] = algo.getBest().getPosition()[i];

  // And test them on our inputs
  odd_numbers.clear();
  even_numbers.clear();
  unclassified_numbers.clear();
  //for(int i = 0 ; i < Problem::nb_digits ; ++i)
  //    printf("%i : %f \n", i, Problem::compute_output(i, params));
  for(size_t i = 0 ; i < p.nb_digits; ++i)
    {
      if(p.compute_output(i, params) >= 0.85)
	odd_numbers.push_back(i);
      else if(p.compute_output(i, params) <= 0.15)
	even_numbers.push_back(i);
      else
	unclassified_numbers.push_back(i);
    }

  printf("I classified %i numbers as EVEN : \n", even_numbers.size());
  p.display_digits(even_numbers);
  printf("I classified %i numbers as ODD : \n", odd_numbers.size());
  p.display_digits(odd_numbers);
  printf("I was not able to classify %i numbers : \n", unclassified_numbers.size());
  p.display_digits(unclassified_numbers);

  delete[] params;
}

