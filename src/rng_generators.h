#ifndef POPOT_RNG_GENERATORS_H
#define POPOT_RNG_GENERATORS_H

#include <cstdlib>
#include <ctime>
#include <iostream>


#include "mersenne.h"

// For random number generations
// There are some discussions in papers of M. Clerc
// e.g. "Randomness matters"

namespace popot
{
  namespace rng
  {

    /*
      JKISS
      See http://www.cs.ucl.ac.uk/staff/D.Jones/GoodPracticeRNG.pdf
    */
   class JKissRNG
    {
      static unsigned int x;
      static unsigned int y;
      static unsigned int z;
      static unsigned int c;
   
    public:
      static int nb_calls;
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
      static constexpr double RNG_RAND_MAX = 4294967296.0;
#else
      static const double RNG_RAND_MAX = 4294967296.0;
#endif

      static void rng_srand()
      {
	rng_srand(time(NULL));
      }
      static void rng_srand(unsigned long seed)
      {
	x = seed | 1;
	y = seed | 2;
	z = seed | 4;
	c = 0;
	nb_calls = 0;
      }
      static void rng_warm_up(void)
      {
	for(int i = 0 ; i < 10000 ; ++i)
	  rng_rand();
      }
      static unsigned int rng_rand(void)
      {
	nb_calls ++;
	unsigned long long t;
	x = 314527869 * x + 1234567;
	y ^= y << 5; y ^= y >> 7; y ^= y << 22;
	t = 4294584393ULL * z + c; c = t >> 32; z = t;
	return x + y + z;
      }
    };

    unsigned int JKissRNG::x;
    unsigned int JKissRNG::y;
    unsigned int JKissRNG::z;
    unsigned int JKissRNG::c;
    int JKissRNG::nb_calls;

    /*

      KISS
      http://www.helsbreth.org/random/rng_kiss.html

      the idea is to use simple, fast, individually promising
      generators to get a composite that will be fast, easy to code
      have a very long period and pass all the tests put to it.
      The three components of KISS are
      x(n)=a*x(n-1)+1 mod 2^32
      y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
      z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
      The y's are a shift register sequence on 32bit binary vectors
      period 2^32-1;
      The z's are a simple multiply-with-carry sequence with period
      2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127
    */

    class KissRNG32
    {
      static unsigned long long kiss_x;
      static unsigned long long kiss_y;
      static unsigned long long kiss_z;
      static unsigned long long kiss_w;
      static unsigned long long kiss_carry;
      static unsigned long long kiss_k;
      static unsigned long long kiss_m;
      
    public:

      static int nb_calls;
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
      static constexpr unsigned long long RNG_RAND_MAX = 4294967295ULL;
#else
      static const unsigned long long RNG_RAND_MAX = 4294967295ULL;
#endif

      static void rng_srand()
      {
	// Random initialization of the seed
	// To really get different seeds, ensure that the calls
	// are separated by at least 1 s.
	rng_srand(time(NULL));
	nb_calls = 0;
      }

      static void rng_srand(unsigned long seed)
      {
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
	nb_calls = 0;
      }

      static void rng_warm_up(void)
      {
	for(int i = 0 ; i < 10000 ; ++i)
	  rng_rand();
      }

      static unsigned long long rng_rand(void)
      {
	kiss_x = kiss_x * 69069 + 1;
	kiss_y ^= kiss_y << 13;
	kiss_y ^= kiss_y >> 17;
	kiss_y ^= kiss_y << 5;
	kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
	kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
	kiss_z = kiss_w;
	kiss_w = kiss_m;
	kiss_carry = kiss_k >> 30;
	nb_calls ++;
	return kiss_x + kiss_y + kiss_w;
      }
    };
    unsigned long long KissRNG32::kiss_x;
    unsigned long long KissRNG32::kiss_y;
    unsigned long long KissRNG32::kiss_z;
    unsigned long long KissRNG32::kiss_w;
    unsigned long long KissRNG32::kiss_carry;
    unsigned long long KissRNG32::kiss_k;
    unsigned long long KissRNG32::kiss_m;
    int KissRNG32::nb_calls;


    /**
     * C Random number generator
     */
    class CRNG
    {     
    public:
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
      static constexpr unsigned long RNG_RAND_MAX = RAND_MAX;
#else
      static const unsigned long RNG_RAND_MAX = RAND_MAX;
#endif
      static int nb_calls;

      static void rng_srand()
      {
	// Random initialization of the seed
	// To really get different seeds, ensure that the calls
	// are separated by at least 1 s.
	CRNG::rng_srand(time(NULL));
      }

      static void rng_srand(unsigned long seed)
      {
	srand(seed);
	nb_calls = 0;
      }

      static void rng_warm_up(void)
      {
	for(int i = 0 ; i < 10000 ; ++i)
	  rng_rand();
      }

      // Generate a random number in [0 ; RAND_MAX [
      static unsigned long rng_rand(void)
      {
	nb_calls++;
	return rand();
      }
    };
    int CRNG::nb_calls;

  }
}

#endif
