/*
    This file is part of popot.
    Copyright (C) 2014  Jeremy Fix, CentraleSupelec

    Author : Jeremy Fix

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact : Jeremy.Fix@centralesupelec.Fr
*/

#ifndef POPOT_NEIGHBORHOOD_H
#define POPOT_NEIGHBORHOOD_H

#include <vector>
#include <algorithm>
#include "exceptions.h"
#include "maths.h"

namespace popot
{
  namespace PSO
  {

    /*
    namespace particle
    {
      template<typename PROBLEM, typename POSITION_INITIALIZER> class BaseParticle;
    }
    */
    namespace neighborhood
    {

      /**
       * Neighborhood
       * @param PARTICLE The particle type you want to use, e.g. swarm::particle::TraditionalParticle
       * \brief A neighborhood contains a vector of particles belonging to the neighborhood as well as a pointer (not a copy!)
       *        to the best particle. The best particle is the first with the lowest fitness 
       */
      template<typename PARTICLE>
      class Neighborhood
      {
      public:
	typedef typename PARTICLE::BestType BestType;
	typedef PARTICLE InNeighborhoodType;

      protected:
	std::vector<InNeighborhoodType *> _particles;
	BestType *_best_particle;

      public:
	Neighborhood(void) : _best_particle(0) 
	{ }


	Neighborhood(const Neighborhood & other)
	{
	  _best_particle = other._best_particle;
	  _particles.resize(other.size());
	  std::copy(other._particles.begin(), other._particles.end(), _particles.begin());

	}

	virtual ~Neighborhood(void){
	  clear();
	}

	void add(InNeighborhoodType * p)
	{
	  _particles.push_back(p);
	}

	unsigned int size() const
	{
	  return _particles.size();
	}

	void clear()
	{
	  _particles.clear();
	}

	std::vector<InNeighborhoodType *>& get() {
	  return _particles;
	}

	InNeighborhoodType * get(unsigned int i)
	{
	  if(i >= 0 && i < size())
	    return _particles[i];
	  else
	    throw popot::Exception::IndexOutOfRange(i, size());
	}

	BestType * findBest(const std::function<int(BestType&, BestType&)>& compare)
	{
	  if(_particles.size() == 0)
	    throw popot::Exception::FindBestFromEmptyNeighborhood();

	  _best_particle = &(_particles[0]->getBestPosition());
	  for(unsigned int i = 1 ; i < _particles.size() ; ++i)
	    {
	      if(compare(_particles[i]->getBestPosition(), *_best_particle) < 0)
		_best_particle = &(_particles[i]->getBestPosition());
	    }

	  return _best_particle;
	}

	BestType* getBest(void)
	{
	  if(_best_particle == 0)
	    throw popot::Exception::BestParticleNotInitialized();
	  return _best_particle;
	}

	void setBest(BestType* best) {
	  _best_particle = best;
	}

	void print(void)
	{
	  std::cout << "Neighborhood hosting " << size() << " particles : " << std::endl;
	  std::cout << "Best particle " << getBest() << " with fitness : " << getBest()->getFitness() << std::endl;
	  for(int i = 0 ; i < size() ; ++i)
	    std::cout << "Particle " << i << " : " << get(i) << std::endl;
	}
      };

      template<typename PARTICLE>
	typename PARTICLE::BestType * findBest(Neighborhood<PARTICLE>* neighborhood,
					       const std::function<int(typename PARTICLE::BestType&, typename PARTICLE::BestType&)>& compare)
	{
	  auto particles = neighborhood->get();
	  if(particles.size() == 0)
	    throw popot::Exception::FindBestFromEmptyNeighborhood();
	  
	  auto best_particle = &(particles[0]->getBestPosition());
	  for(unsigned int i = 1 ; i < particles.size() ; ++i)
	    {
	      if(compare(particles[i]->getBestPosition(), *best_particle) < 0)
		best_particle = &(particles[i]->getBestPosition());
	    }
	  neighborhood->setBest(best_particle);

	  return best_particle;
	}

    } // namespace neighborhood
  } // namespace PSO
} // namespace popot

#endif // POPOT_NEIGHBORHOOD_H
