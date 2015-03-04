// In this example, we show how to make use of the random number generators provided by the library

#include "rng_generators.h"

typedef popot::rng::JKissRNG RNG_GENERATOR;
#define FILENAME "jkissrng.data"

/*
typedef popot::rng::CRNG RNG_GENERATOR;
#define FILENAME "crng.data"
*/

#include "popot.h"

int main(int argc, char * argv[])
{
  std::ofstream outfile(FILENAME);

  // Let's initialize the seed
  RNG_GENERATOR::rng_srand();

  // And generate some samples
   for(int i = 0 ; i < 100 ; ++i)
     outfile << i << " " << popot::math::uniform_random(0,1) << std::endl;
   outfile.close();

}
