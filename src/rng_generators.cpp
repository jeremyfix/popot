#include "rng_generators.h"

unsigned int popot::rng::JKissRNG::x;
unsigned int popot::rng::JKissRNG::y;
unsigned int popot::rng::JKissRNG::z;
unsigned int popot::rng::JKissRNG::c;
int popot::rng::JKissRNG::nb_calls;


void popot::rng::JKissRNG::rng_srand() {
  popot::rng::JKissRNG::rng_srand(time(NULL));
}

void popot::rng::JKissRNG::rng_srand(unsigned long seed)  {
  x = seed | 1;
  y = seed | 2;
  z = seed | 4;
  c = 0;
  nb_calls = 0;
}

void popot::rng::JKissRNG::rng_warm_up(void) {
  for(int i = 0 ; i < 10000 ; ++i)
    rng_rand();
}

unsigned int popot::rng::JKissRNG::rng_rand(void) {
  ++nb_calls;
  unsigned long long t;
  x = 314527869 * x + 1234567;
  y ^= y << 5; y ^= y >> 7; y ^= y << 22;
  t = 4294584393ULL * z + c; c = t >> 32; z = t;
  return x + y + z;
}

unsigned long long popot::rng::KissRNG32::kiss_x;
unsigned long long popot::rng::KissRNG32::kiss_y;
unsigned long long popot::rng::KissRNG32::kiss_z;
unsigned long long popot::rng::KissRNG32::kiss_w;
unsigned long long popot::rng::KissRNG32::kiss_carry;
unsigned long long popot::rng::KissRNG32::kiss_k;
unsigned long long popot::rng::KissRNG32::kiss_m;
int popot::rng::KissRNG32::nb_calls;

void popot::rng::KissRNG32::rng_srand()  {
  // Random initialization of the seed
  // To really get different seeds, ensure that the calls
  // are separated by at least 1 s.
  popot::rng::KissRNG32::rng_srand(time(NULL));
}

void popot::rng::KissRNG32::rng_srand(unsigned long seed) {
  kiss_x = seed | 1;
  kiss_y = seed | 2;
  kiss_z = seed | 4;
  kiss_w = seed | 8;
  kiss_carry = 0;
  nb_calls = 0;
}

void popot::rng::KissRNG32::rng_warm_up(void) {
  for(int i = 0 ; i < 10000 ; ++i)
    rng_rand();
}

unsigned long long popot::rng::KissRNG32::rng_rand(void) {
  kiss_x = kiss_x * 69069 + 1;
  kiss_y ^= kiss_y << 13;
  kiss_y ^= kiss_y >> 17;
  kiss_y ^= kiss_y << 5;
  kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
  kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
  kiss_z = kiss_w;
  kiss_w = kiss_m;
  kiss_carry = kiss_k >> 30;
  ++nb_calls;
  return kiss_x + kiss_y + kiss_w;
}

int popot::rng::CRNG::nb_calls;

void popot::rng::CRNG::rng_srand() {
  // Random initialization of the seed
  // To really get different seeds, ensure that the calls
  // are separated by at least 1 s.
  popot::rng::CRNG::rng_srand(time(NULL));
}

void popot::rng::CRNG::rng_srand(unsigned long seed) {
  srand(seed);
  nb_calls = 0;
}

void popot::rng::CRNG::rng_warm_up(void) {
  for(int i = 0 ; i < 10000 ; ++i)
    rng_rand();
}


unsigned long popot::rng::CRNG::rng_rand(void) {
  ++nb_calls;
  return rand();
}
