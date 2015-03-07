/*
  This file is part of popot.

  Copyright (C) 2014, Jeremy Fix, CentraleSupelec

  Author : Jeremy Fix

  popot is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  popot is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with popot.  If not, see <http://www.gnu.org/licenses/>.

  Contact : Jeremy.Fix@centralesupelec.fr
*/

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
      static constexpr double RNG_RAND_MAX = 4294967296.0;
      
      static void rng_srand();
      static void rng_srand(unsigned long seed);

      static void rng_warm_up(void);
      static unsigned int rng_rand(void);
    };

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
      static constexpr unsigned long long RNG_RAND_MAX = 4294967295ULL;

      static void rng_srand();
      static void rng_srand(unsigned long seed);
      static void rng_warm_up(void);
      static unsigned long long rng_rand(void);
    };


    /**
     * C Random number generator
     */
    class CRNG
    {     
    public:
      static constexpr unsigned long RNG_RAND_MAX = RAND_MAX;
      static int nb_calls;

      static void rng_srand();
      static void rng_srand(unsigned long seed);
      static void rng_warm_up(void);
      // Generate a random number in [0 ; RAND_MAX [
      static unsigned long rng_rand(void);
    };

  }
}

#endif
