#ifndef POPOT_TOPOLOGY_H
#define POPOT_TOPOLOGY_H

#include <vector>
#include <map>
#include "neighborhood.h"

// Define a connection matrix 
// and make a generic method which takes this connection matrix;
// and fills in the neighborhoods

namespace popot
{
  namespace PSO
  {
    namespace topology
    {

      /**
       *
       * Basic topology methods
       */
      
      template<typename PARTICLE>
      void connectParticles(std::vector<PARTICLE>& particles,
			    std::map< size_t, std::vector<size_t> > &neighbordhood_membership,
			    bool * who_informs_whom)
      {
	size_t size = particles.size();

	// Just clean up the neighborhood memberships
	for(size_t i = 0 ; i < size ; ++i)
	  neighbordhood_membership[i].clear();



	for(size_t j = 0 ; j < size ; ++j)
	  {
	    particles[j].getNeighborhood().clear();

	    // We browse column j of who_informs_whom
	    // if 1.0, then we add particle i to the neighborhood of particle j
	    for(size_t i = 0 ; i < size ; ++i)
	      // If particle i informs particle j
	      if(who_informs_whom[i*size + j])
		{
		  // Particle i belongs to the neighborhood of particle j
		  neighbordhood_membership[i].push_back(j);
		  // And we add particle i to the neighborhood of particle j
		  particles[j].getNeighborhood().add(&(particles[i]));
		}
	  }
      }

      template<typename PARTICLE>
      void getNeighborhoodList(std::vector<PARTICLE>& particles,
			       std::vector< typename PARTICLE::NeighborhoodType *> &neighborhoods)
      {
	neighborhoods.clear();
	for(size_t i = 0 ; i < particles.size() ; ++i)
	  neighborhoods.push_back(&(particles[i].getNeighborhood()));
      }

      /**
       * Full toplogy
       * @short Each of the N particles receives information from the N others
       */
      template<typename PARTICLE, bool SELF>
      void full_fillNeighborhoods(std::vector<PARTICLE> &particles,
				  std::vector< typename PARTICLE::NeighborhoodType *> &neighborhoods,
				  std::map< size_t, std::vector<size_t> > &neighbordhood_membership)
      {
	size_t size = particles.size();

	// Generate the connection matrix
	// Column j indicates which particle informs particle j
	// therefore, line i indicates which particles the particle i informs
	bool *who_informs_whom = new bool[size * size];

	memset(who_informs_whom, true, size*size*sizeof(bool));

	if(!SELF)
	  for(size_t i = 0 ; i < size ; ++i)
	    who_informs_whom[i*size + i] = false;

	// Given the connectivity matrix, we now connect the particles
	getNeighborhoodList(particles, neighborhoods);
	connectParticles(particles, neighbordhood_membership, who_informs_whom);
	  
	delete[] who_informs_whom;
      }

      /**
       * Ring toplogy
       * @short A ring topology connects the particles on a ring plus a connection to itself
       */
      template<typename PARTICLE, bool SELF>
      void ring_fillNeighborhoods(std::vector<PARTICLE>& particles,
				  std::vector< typename PARTICLE::NeighborhoodType *> &neighborhoods,
				  std::map< size_t, std::vector<size_t> > &neighbordhood_membership)
      	{
	  size_t size = particles.size();
	  
      	  // Generate the connection matrix
      	  // Column j indicates which particle informs particle j
      	  // therefore, line i indicates which particles the particle i informs
      	  bool *who_informs_whom = new bool[size * size];

      	  // Set the connection matrix to False everywhere
	  memset(who_informs_whom, false, size*size*sizeof(bool));

       	  size_t i_neigh;
      	  for(size_t i = 0 ; i < size ; ++i)
      	    {
      	      // Particle i is informed by the particle on the left and the particle on its right
      	      if(i == size - 1)
		i_neigh = 0;
	      else
		i_neigh = i + 1;
      	      who_informs_whom[i_neigh*size + i] = true;

      	      if(i == 0)
      		i_neigh = size-1;
	      else
		i_neigh = i - 1;
      	      who_informs_whom[i_neigh*size + i] = true;

      	      if(SELF)
      		who_informs_whom[i*size + i] = true;
      	    }

      	  // Given the connectivity matrix, we now connect the particles
      	  getNeighborhoodList<PARTICLE>(particles, neighborhoods);
      	  connectParticles<PARTICLE>(particles, neighbordhood_membership, who_informs_whom);

      	  delete[] who_informs_whom;
      	}
      	

      /**
       * Von Neuman topology
       * @short The von Neuman topology connects the particles on a 2D toric grid
       *        with the nearest four neighbours plus itself
       */
      template<typename PARTICLE, bool SELF>
      void vonNeuman_fillNeighborhoods(std::vector<PARTICLE>& particles,
				       std::vector< typename PARTICLE::NeighborhoodType *> &neighborhoods,
				       std::map< size_t, std::vector<size_t> > &neighbordhood_membership)
      {
	size_t size = particles.size();
	size_t width = size_t(sqrt(size));
	size_t height= size_t(size / width);

	if(height * width != size)
	  {
	    std::cout << "WARNING : I was only able to set up a grid of " << height << " x " << width << ", not including all the " << size << " particles " << std::endl;
	    std::cout << " this will fail ;) " << std::endl;
	    
	  }

	// Generate the connection matrix
	bool *who_informs_whom = new bool[size * size];

	// Column j indicates which particle informs particle j
	// therefore, line i indicates which particles the particle i informs

	// Set the connection matrix to 0 everywhere
	for(size_t i = 0 ; i < size * size ; ++i)
	  who_informs_whom[i] = false;

	size_t i_neigh, j_neigh;
	size_t index_part, index_neigh;
	for(size_t i = 0 ; i < height ; ++i)
	  {
	    for(size_t j = 0 ; j < width ; ++j)
	      {

		index_part = i*width + j;

		// Add the 4 closest neighboors that are informing particle index_part
		if(i == 0)
		  i_neigh = height - 1;
		else
		  i_neigh = i - 1;
		j_neigh = j;

		index_neigh = i_neigh * width + j_neigh;
		who_informs_whom[index_neigh*size + index_part] = true;

		if(i == height - 1)
		  i_neigh = 0;
		else
		  i_neigh = i + 1;
		j_neigh = j;

		index_neigh = i_neigh * width + j_neigh;
		who_informs_whom[index_neigh*size + index_part] = true;

		i_neigh = i;
		if(j == 0)
		  j_neigh = width - 1;
		else
		  j_neigh = j - 1;

		index_neigh = i_neigh * width + j_neigh;
		who_informs_whom[index_neigh*size + index_part] = true;

		i_neigh = i;
		if(j_neigh == width - 1)
		  j_neigh = 0;
		else
		  j_neigh = j + 1;

		index_neigh = i_neigh * width + j_neigh;
		who_informs_whom[index_neigh*size + index_part] = true;

		if(SELF)
		  who_informs_whom[i*size + i] = true;
	      }
	  }
	
	// Given the connectivity matrix, we now connect the particles
	getNeighborhoodList<PARTICLE>(particles, neighborhoods);
	connectParticles<PARTICLE>(particles, neighbordhood_membership, who_informs_whom);

	delete[] who_informs_whom;
      }



      /**
       * Random topology
       * @short Some probabilistic connectivity so that most of the particles will have K informants
       */
      template<typename PARTICLE, int K, bool SELF>
      void randomInformants_fillNeighborhoods(std::vector<PARTICLE>& particles,
					      std::vector< typename PARTICLE::NeighborhoodType *>& neighborhoods,
					      std::map< size_t, std::vector<size_t> >& neighbordhood_membership)
      	{

	  size_t size = particles.size();

      	  // Generate the connection matrix
      	  // Column j indicates which particle informs particle j
      	  // therefore, line i indicates which particles the particle i informs
      	  bool *who_informs_whom = new bool[size * size];

      	  // Set the connection matrix to 0 everywhere
      	  for(size_t i = 0 ; i < size*size ; ++i)
      	    who_informs_whom[i] = false;

      	  // A particle informs itself
      	  if(SELF)
      	    for(size_t i = 0 ; i < size ; ++i)
      	      who_informs_whom[i*size + i] = true;

      	  for(size_t i = 0 ; i < size ; ++i)
      	    {
      	      // For line i, we draw uniformely with replacement K-1 particles this particle will inform
      	      for(size_t k = 0 ; k < K-1 ; ++k)
      		who_informs_whom[i*size + (size_t)popot::math::uniform_integer(0, size-1)] = true;
      	    }

      	  // Given the connectivity matrix, we now connect the particles
	  getNeighborhoodList<PARTICLE>(particles, neighborhoods);
	  connectParticles<PARTICLE>(particles, neighbordhood_membership, who_informs_whom);

      	  delete[] who_informs_whom;
      	}

      /**
       * Adaptive random topology
       * @short Some probabilistic connectivity so that most of the particles will have K informants
       */
      template<typename PARTICLE, int K, bool SELF>
      void adaptiveRandom_fillNeighborhoods(std::vector<PARTICLE>& particles,
      			       std::vector< typename PARTICLE::NeighborhoodType *> &neighborhoods,
			     std::map< size_t, std::vector<size_t> > &neighbordhood_membership)
      	{
	  size_t size = particles.size();

      	  double p = 1.0 - pow(1.0 - 1.0/size,K);

      	  // Generate the connection matrix
      	  // Column j indicates which particle informs particle j
      	  // therefore, line i indicates which particles the particle i informs
      	  bool *who_informs_whom = new bool[size * size];
	  memset(who_informs_whom, false, size * size * sizeof(bool));

      	  // A particle informs itself
      	  if(SELF)
      	    for(size_t i = 0 ; i < size ; ++i)
      	      who_informs_whom[i*size + i] = true;

      	  for(size_t i = 0 ; i < size ; ++i)
      	    {
      	      for(size_t k = 0 ; k < size ; ++k)
      		if(k != i)
      		  {
      		    if(popot::math::uniform_random(0.0,1.0) < p)
      		      who_informs_whom[k*size + i] = true; // k informs i
      		    else
      		      who_informs_whom[k*size + i] = false;
      		  }
      	    }

      	  // Given the connectivity matrix, we now connect the particles
	  getNeighborhoodList<PARTICLE>(particles, neighborhoods);
	  connectParticles<PARTICLE>(particles, neighbordhood_membership, who_informs_whom);

      	  delete[] who_informs_whom;
      	}



    } // topology
  } // PSO
} // popot

#endif // POPOT_TOPOLOGY_H
