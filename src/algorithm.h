#ifndef POPOT_ALGORITHM_H
#define POPOT_ALGORITHM_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "neighborhood.h"
#include "individuals.h"
#include "topology.h"
#include "confinements.h"
#include <fstream>


// TODO
// - Rendre encore plus générique les manipulations de vectors : tout vectoriel!
//   et utiliser des itérateurs plutôt que [i]
// - ajouter les options -ansi -pedantic : catastrophe, ça pète sur les auto ... en tout cas sur fc16

namespace popot
{
  namespace PSO
  {
    namespace algorithm
    {
      typedef enum
	{
	  SYNCHRONOUS_EVALUATION,
	  ASYNCHRONOUS_EVALUATION,
	  ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION
	} EvaluationMode;

      /**
       * Standard particle swarm optimization
       * @param PARTICLE The particle type you want to use, e.g. swarm::particle::TraditionalParticle
       * @param TOPOLOGY Defines the topology, e.g. swarm::topology::VonNeuman
       * @param STOP_CRITERIA Defines the condition to stop the evolution of the swarm
       */

      template< typename LBOUND_FUNC, typename UBOUND_FUNC, 
		typename STOP_CRITERIA, typename COST_FUNCTION, 
		typename TOPOLOGY, typename UPDATE_POSITION_RULE, 
		typename UPDATE_VELOCITY_RULE, typename UPDATE_BEST_POSITION_RULE, 
		typename CONFINE_FUNCTION, typename INIT_POSITION_FUNCTION, 
		typename INIT_VELOCITY_FUNCTION, typename PARTICLE>
      class Base
      {
      public:
	typedef typename PARTICLE::BestType BestType;
	typedef PARTICLE  ParticleType;
	typedef typename PARTICLE::NeighborhoodType NeighborhoodType;

      private:
	std::vector<PARTICLE> particles;
	std::vector<NeighborhoodType * > neighborhoods;
	std::map< size_t , std::vector<size_t> > neighborhood_membership;

	size_t *particles_indexes;

	BestType _best_particle;

	size_t _dimension;
	size_t _swarm_size;

	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;
	const STOP_CRITERIA& _stop_criteria;
	const COST_FUNCTION& _cost_function;
	const TOPOLOGY _topology;
	const UPDATE_POSITION_RULE _update_position_rule;
	const UPDATE_VELOCITY_RULE _update_velocity_rule;
	const UPDATE_BEST_POSITION_RULE _update_best_position_rule;

	const CONFINE_FUNCTION _confine_function;
	const INIT_POSITION_FUNCTION _init_position_function;
	const INIT_VELOCITY_FUNCTION _init_velocity_function;

	const bool _reevaluate_best_before_updating;

	EvaluationMode _evaluation_mode;
      public:
	int epoch;
	int nb_new_neigh;

      public:
	Base(size_t swarm_size,
	     size_t dimension,
	     const LBOUND_FUNC& lbound,
	     const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop,
	     const COST_FUNCTION& cost_function,    
	     const TOPOLOGY topology,
	     const UPDATE_POSITION_RULE update_position_rule,
	     const UPDATE_VELOCITY_RULE update_velocity_rule,
	     const UPDATE_BEST_POSITION_RULE update_best_position_rule,
	     const CONFINE_FUNCTION confine_function,
	     const INIT_POSITION_FUNCTION init_position_function,
	     const INIT_VELOCITY_FUNCTION init_velocity_function,
	     const PARTICLE& p,
	     EvaluationMode evaluation_mode,
	     bool reevaluate_best_before_updating) 
	: particles_indexes(0),
	  _best_particle(dimension),
	  _dimension(dimension),
	  _swarm_size(swarm_size),
	  _lbound(lbound),
	  _ubound(ubound),
	  _stop_criteria(stop),
	  _cost_function(cost_function),
	  _topology(topology),
	  _update_position_rule(update_position_rule),
	  _update_velocity_rule(update_velocity_rule),
	  _update_best_position_rule(update_best_position_rule),
	  _confine_function(confine_function),
	  _init_position_function(init_position_function),
	  _init_velocity_function(init_velocity_function),
	  _reevaluate_best_before_updating(reevaluate_best_before_updating),
	  _evaluation_mode(evaluation_mode),
	  epoch(0),
	  nb_new_neigh(0)
	{
	  init();
	}

	virtual ~Base(void)
	{
	  delete[] particles_indexes;
	}

	void init(void)
	{
	
	  epoch = 0;
	  nb_new_neigh = 0;
	  //
	  if(particles_indexes != 0)
	    delete[] particles_indexes;
	  particles_indexes = new size_t[_swarm_size];

	  // Declare our particles
	  particles.clear();
	  for(size_t i = 0 ; i < _swarm_size ; ++i)
	    {
	      // Create the particle
	      particles.push_back(PARTICLE(_dimension));

	      // Initialize the position, best position and velocity
	      _init_position_function(particles[i].getPosition(), _lbound, _ubound);
	      _init_velocity_function(particles[i].getPosition(), particles[i].getVelocity(), _lbound, _ubound);
				       
	      // Evaluate the initial fitnesses
	      particles[i].evaluateFitness(_cost_function);
		
	      // And initialize the best position
	      particles[i].initBestPosition();

	      // Set up the array of indices, used to random shuffle 
	      particles_indexes[i] = i;
	    }
	    
	  // We now form groups of particles depending on the topology
	  _topology(particles, neighborhoods, neighborhood_membership);
	  nb_new_neigh ++;

	  // We now browse all the neighborhoods and find the best particles
	  // within each of them
	  // and initialize the best best particle of the whole swarm
	  _best_particle = *(neighborhoods[0]->findBest());
	  for(size_t i = 1 ; i < neighborhoods.size() ; ++i)
	    if(neighborhoods[i]->findBest()->compare(_best_particle) < 0)
	      _best_particle = *(neighborhoods[i]->getBest());
	}


	/**
	 * Returns the size of the swarm
	 */
	size_t getSize(void)
	{
	  return _swarm_size;
	}

	/**
	 * One step of the swarm
	 */
	double step(void)
	{
	  // For each particle,
	  // - Update the particle's velocity
	  // - Enforce velocity boundaries
	  // - Move the particle to its new position
	  // - Enforce position boundaries 
	  // - Update the particle's best position
	  // - Update the swarm's best position

	  // if(VERBOSE_BENCH)
	  //   std::cout << "Random before looping : " << RNG_GENERATOR::nb_calls << " calls " << std::endl;

	  // If we use an asynchronous update, we first shuffle the particles
	  if(_evaluation_mode == ASYNCHRONOUS_EVALUATION || _evaluation_mode == ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION)
	    {
	      // If comparing to SPSO, the following should be commented
	      if(_evaluation_mode == ASYNCHRONOUS_EVALUATION)
	      	popot::math::random_shuffle_indexes(particles_indexes, _swarm_size);

	      size_t particle_index;
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
	   	{
		  //std::cout << "updating particle " << i << std::endl;
	   	  particle_index = particles_indexes[i];

	   	  // First find the best informant within the neighborhood
	   	  particles[particle_index].getNeighborhood().findBest();

	   	  // Update the velocities and position
		  _update_velocity_rule(particles[particle_index]);
	   	  _update_position_rule(particles[particle_index]);

	   	  // Confine the positions and velocities if required
		  _confine_function(particles[particle_index].getPosition(), particles[particle_index].getVelocity(), _lbound, _ubound);
		  //std::cout << "fitness " << i << std::endl;
	   	  // Compute the fitness of the new position
	   	  particles[particle_index].evaluateFitness(_cost_function);
		  //std::cout << "ok " << i << std::endl;
	   	  // And see if we update the personal best
		  if(_reevaluate_best_before_updating)
		    particles[i].getBestPosition().evaluateFitness(_cost_function);
		  //std::cout << "update best " << i << std::endl;
		  _update_best_position_rule(particles[particle_index]);

		  //std::cout << "looping " << i << std::endl;
	   	  //if(VERBOSE_BENCH)
	   	  //  std::cout << std::endl;

	   	}
	    }
	  else // Synchronous evaluation
	    {
	      // In synchronous mode
	      // we first update all the current positions and evaluate their fitness
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
	  	{
	   	  // Update the velocities and position
		  _update_velocity_rule(particles[i]);
	   	  _update_position_rule(particles[i]);

	   	  // Confine the positions and velocities if required
		  _confine_function(particles[i].getPosition(), particles[i].getVelocity(), _lbound, _ubound);
	   	  // Compute the fitness of the new position
	   	  particles[i].evaluateFitness(_cost_function);
	  	}
	      // before changing the personal bests


	      // The personal best position is updated after all the positions of the particles
	      // have been updated because the neighborhoods hold a pointer to a personal best
	      // As this best in the neighborhood is used by the particles when they
	      // update their velocity, we must ensure that, for a synchronous evaluation
	      // updating the best is done after all the positions updates
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
		{
		  if(_reevaluate_best_before_updating)
		    particles[i].getBestPosition().evaluateFitness(_cost_function);
		  _update_best_position_rule(particles[i]);
		}

	      // We now update the best particle of all the neighborhoods
	      for(size_t i = 0 ; i < neighborhoods.size() ; ++i)
	  	neighborhoods[i]->findBest();
	    }

	  // Update the best particle the whole swarm ever had
	  double old_fitness = _best_particle.getFitness();

	  _best_particle = *(neighborhoods[0]->findBest());
	  for(size_t i = 1 ; i < neighborhoods.size() ; ++i)
	    if(neighborhoods[i]->findBest()->compare(_best_particle) < 0)
	      _best_particle = *(neighborhoods[i]->getBest());

	  if(_best_particle.getFitness() >= old_fitness)
	    {
	      // We consider that there is no improvement in the best particle
	      // and ask the topology if it wants to regenerate its topology
	      _topology(particles, neighborhoods, neighborhood_membership);
	      nb_new_neigh ++;
	    }
	  
	  epoch++;
	  return _best_particle.getFitness() - old_fitness;
	}

	/**
	 * Current epoch
	 */
	size_t getEpoch(void) const
	{
	  return epoch;
	}

	/**
	 * Run until stop criteria is met. 
	 * @param verbose set it to 1 to print the best fitness at each iteration
	 */
	void run(int verbose=0)
	{
	  while(!_stop_criteria(getBest().getFitness(), epoch))
	    {
	      step();
	      if(verbose) std::cout << '\r' << std::setw(6) << std::setfill('0') << epoch << " " << getBest().getFitness() << std::setw(5) << std::setfill(' ') << ' ' << std::flush;
	    }
	  if(verbose) std::cout << std::endl;
	}

	/**
	 * Returns the current best particle of the whole swarm (best of the personal best)
	 */
	BestType& getBest(void)
	{
	  return _best_particle;
	}

	/**
	 * Returns the fitness of the best individual (Required for benchmarking)
	 */ 
	double getBestFitness() const
	{
	  return _best_particle.getFitness();
	}

	/**
	 * Returns a reference to the vector of particles
	 */
	std::vector<PARTICLE>& getParticles(void)
	{
	  return particles;
	}

	/**
	 * Print in the console the status of the swarm
	 * @param mode 0 displays all the particles, 1 display only the best
	 */
	void print(int mode=0)
	{
	  if(_evaluation_mode == SYNCHRONOUS_EVALUATION)
	    printf("Using a synchronous evaluation mode \n");
	  else if(_evaluation_mode == ASYNCHRONOUS_EVALUATION)
	    printf("Using an asynchronous evaluation mode \n");
	  else if(_evaluation_mode == ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION)
	    printf("Using an asynchronous evaluation mode with random shuffling the particles\n");
	  else
	    printf("WARNING : Unrecognized evaluation mode \n");

	  switch(mode)
	    {
	    case 0:
	      {
		for(size_t i = 0 ; i < _swarm_size ; ++i)
		  {
		    std::cout << "Particle " << std::setw(3) << std::setfill(' ') << i
			      << " : " << "\n" << particles[i] << std::endl;
		  }
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle << std::endl;
	      }
	      break;
	    case 1:
	      {
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle << std::endl;
	      }
	      break;
	    case 2:
	      {
		for(size_t i = 0 ; i < _swarm_size ; ++i)
		  {
		    std::cout << "Particle " << std::setw(3) << std::setfill(' ') << i
			      << " : " << particles[i].getFitness() << " ; Best : " << particles[i].getBestPosition().getFitness() << std::endl;
		  }
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle.getFitness() << std::endl;
	      }
	      break;
	    default:
	      break;
	    }
	}

	
	std::map< size_t , std::vector<size_t> >& getNeighborhoodMembership(void)
	{
	  return neighborhood_membership;
	}

	/**
	 * Generate a graph of the connections between the particles using DOT
	 */
	void generateGraph(std::string filename)
	{
	  printf("----------- \n Generating a graph of the connectivity \n");
	  size_t nb_particles = _swarm_size;

	  // From this connectivity matrix, we generate the file
	  // to be used with the DOT command
	  std::ofstream outfile(filename.c_str());
	  outfile << "digraph Connections {" << std::endl;
	  outfile << "node [shape=circle,fixedsize=true,width=0.6];  ";
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      outfile << i << ";";
	    }
	  outfile << std::endl;

	  // Generate a connectivity matrix
	  int *connectivity_matrix = new int[nb_particles*nb_particles];
	  for(size_t i = 0 ; i < nb_particles*nb_particles ; ++i)
	    connectivity_matrix[i] = 0;

	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    for(size_t j = 0 ; j < neighborhood_membership[i].size(); ++j)
	      connectivity_matrix[i*nb_particles+neighborhood_membership[i][j]] += 1;

	  // Put in the connections
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      for(size_t j = 0 ; j < i ; ++j)
		{
		  if(connectivity_matrix[i*nb_particles+j] && !connectivity_matrix[j*nb_particles+i])
		    outfile << j << "->" << i << ";" << std::endl;
		  else if(connectivity_matrix[i*nb_particles+j] && connectivity_matrix[j*nb_particles+i])
		    outfile << j << "->" << i << " [dir=both];" << std::endl;
		}
	      if(connectivity_matrix[i*nb_particles+i])
		outfile << i << "->" << i << ";" << std::endl;
	    }

	  // The title and the end of the graphic
	  outfile << "overlap=false" << std::endl;
	  outfile << "label=\"Swarm's connections\"" << std::endl;
	  outfile << "fontsize=12;" << std::endl;
	  outfile << "}" << std::endl;
	  outfile.close();

	  printf("You can generate the graph by calling e.g.: \n");
	  printf(" neato -Tpng %s > graph.png\n", filename.c_str());
	  printf("Other export formats can be used with dot/neato : svg, eps, dia, xfig, ...\n");

	  // Compute the distribution of the informants
	  // We first compute the number of informants for each particle
	  // The informants of particle i are at line i
	  size_t * nb_informants = new size_t[nb_particles];
	  size_t max_nb_informants = 0;
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      nb_informants[i] = 0;
	      for(size_t j = 0 ; j < nb_particles ; ++j)
		nb_informants[i] += connectivity_matrix[i*nb_particles+j];
	      if(nb_informants[i] > max_nb_informants)
		max_nb_informants = nb_informants[i];
	    }
	  printf("\n");
	  printf("Max number of informants : %i \n", max_nb_informants);

	  outfile.open((filename + ".conn").c_str());
	  size_t nb_part = 0;
	  for(size_t i = 0 ; i <= max_nb_informants ; ++i)
	    {
	      // How many particles have nb_informants <= i ?
	      nb_part = 0;
	      for(size_t j = 0 ; j < nb_particles ; ++j)
		if(nb_informants[j] == i)
		  nb_part ++;
	      outfile << i << "\t" << nb_part/double(nb_particles) << std::endl;
	    }
	  outfile.close();
	  printf("-------------\n");

	  delete[] connectivity_matrix;
	  delete[] nb_informants;
	}

	void save(const char* filename) {
	  std::cout << "Saving the particles in " << filename << std::endl;
	  std::ofstream outfile(filename);

	  // Iterate over the particles and ask them to dump
	  // their content into outfile
	  for(size_t i = 0 ; i < _swarm_size ; ++i)
	      particles[i].save(outfile);

	  outfile.close();
	  	  
	}

	void load(const char* filename) {
	  std::cout << "Loading the particles from " << filename << std::endl;
	  std::ifstream infile(filename);
	  for(size_t i = 0 ; i < _swarm_size; ++i) {
	    particles[i].load(infile);
	    particles[i].evaluateFitness(_cost_function);
	    particles[i].initBestPosition();

	  }
	  infile.close();

	  // We now browse all the neighborhoods and find the best particles
	  // within each of them
	  // and initialize the best best particle of the whole swarm
	  _best_particle = *(neighborhoods[0]->findBest());
	  for(size_t i = 1 ; i < neighborhoods.size() ; ++i)
	    if(neighborhoods[i]->findBest()->compare(_best_particle) < 0)
	      _best_particle = *(neighborhoods[i]->getBest());


	}
      };

      template< typename LBOUND_FUNC, typename UBOUND_FUNC,
		typename STOP_CRITERIA, typename COST_FUNCTION,
		typename TOPOLOGY, typename UPDATE_POSITION_RULE, 
		typename UPDATE_VELOCITY_RULE,typename UPDATE_BEST_POSITION_RULE, 
		typename CONFINE_FUNCTION, typename INIT_POSITION_FUNCTION,
		typename INIT_VELOCITY_FUNCTION, typename PARTICLE>
      Base<LBOUND_FUNC, UBOUND_FUNC, 
	   STOP_CRITERIA, COST_FUNCTION, 
	   TOPOLOGY, UPDATE_POSITION_RULE, 
	   UPDATE_VELOCITY_RULE, UPDATE_BEST_POSITION_RULE, 
	   CONFINE_FUNCTION, INIT_POSITION_FUNCTION, 
	   INIT_VELOCITY_FUNCTION, PARTICLE>
      base(size_t swarm_size,
	   size_t dimension,
	   const LBOUND_FUNC& lbound,
	   const UBOUND_FUNC& ubound,
	   const STOP_CRITERIA& stop,
	   const COST_FUNCTION& cost_function,
	   const TOPOLOGY& topology,
	   const UPDATE_POSITION_RULE& update_position_rule,
	   const UPDATE_VELOCITY_RULE& update_velocity_rule,
	   const UPDATE_BEST_POSITION_RULE& update_best_position_rule,
	   const CONFINE_FUNCTION& confine_function,
	   const INIT_POSITION_FUNCTION& init_position_function,
	   const INIT_VELOCITY_FUNCTION& init_velocity_function,
	   const PARTICLE& p,
	   EvaluationMode evaluation_mode,
	   bool reevaluate_best_before_updating)
      {
	// the parameters above can be provided by reference
	// as we anyway copy (some of) them in the constructor
	
	return Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
		    TOPOLOGY, UPDATE_POSITION_RULE, UPDATE_VELOCITY_RULE, UPDATE_BEST_POSITION_RULE, CONFINE_FUNCTION, INIT_POSITION_FUNCTION, INIT_VELOCITY_FUNCTION, PARTICLE>
	  (swarm_size, dimension, lbound, ubound, stop, cost_function, topology, 
	   update_position_rule, update_velocity_rule, update_best_position_rule, confine_function, init_position_function, init_velocity_function, p, evaluation_mode, reevaluate_best_before_updating);
      }


    } // namespace algorithm
  } // namespace PSO

  namespace ABC
  {
    namespace algorithm
    {
      template<typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      class Base
      {
	typedef popot::ABC::individuals::FoodSource FoodSourceType;

      private:
	size_t _epoch;
	size_t _CS;
	size_t _dimension;
	size_t _limitForScout;
	size_t _nb_employed;
	size_t _nb_onlookers;
	double * _probabilities;

	FoodSourceType * _foodSources;
	FoodSourceType _bestSource;

	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;
	const STOP_CRITERIA& _stop_criteria;
	const COST_FUNCTION& _cost_function;

	
      private:
	void findBestSource(void)
	{
	  double bestf = _bestSource.getFitness();
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    if(bestf < _foodSources[i].getFitness())
	      {
		bestf = _foodSources[i].getFitness();
		_bestSource = _foodSources[i];
	      }
	}

	void employedPhase(void)
	{
	  // In the employed bees phase, we generate new solutions
	  // around each nectar source
	  //int change_dim;
	  size_t other_source;
	  //double phi;
	  FoodSourceType new_source;
	  //double new_param_value;
	  double sum_fitnesses = 0;
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    {
	      // Randomly select another source, different from the current source
	      other_source = (size_t) popot::math::uniform_random(0, _nb_employed);
	      while(other_source == i)
		other_source = (size_t) popot::math::uniform_random(0, _nb_employed);

	      _foodSources[i].combine(_foodSources[other_source], _lbound, _ubound, _cost_function);

	      _probabilities[i] = (_foodSources[i]).getFitness();
	      sum_fitnesses += _probabilities[i];
	    }

	  // At the end, we normalize the probabilities for each source
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _probabilities[i] /= sum_fitnesses;
	}

	void onlookerPhase(void)
	{
	  // Onlooker phase
	  // Each onlooker bee selects a food source
	  // based on its probability (reflecting its relative fitness)
	  size_t selected_source=0;
	  size_t other_source=0;
	  for(size_t i = 0 ; i < _nb_onlookers ; ++i)
	    {
	      // Select a source based on its fitness
	      selected_source = popot::math::random_from_array(_probabilities);

	      // Randomly select another source, different from the current source
	      other_source = (size_t) popot::math::uniform_random(0, _nb_employed);
	      while(other_source == i)
		other_source = (size_t) popot::math::uniform_random(0, _nb_employed);

	      _foodSources[selected_source].combine(_foodSources[other_source], _lbound, _ubound, _cost_function);
	      
	    }
	}

	void scoutPhase(void)
	{
	  // Scout phase
	  // We browse all the food sources
	  // If a source has a counter higher than the limit
	  // we reset it
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    if(_foodSources[i].getCounter() >= _limitForScout)
	      _foodSources[i].init(_lbound, _ubound, _cost_function);
	}

      public:

	Base(int colony_size, int dimension, const LBOUND_FUNC &lbound, const UBOUND_FUNC &ubound, const STOP_CRITERIA &stop_criteria, const COST_FUNCTION &cost_function)
	  : _epoch(0), 
	    _CS(colony_size), 
	    _dimension(dimension), 
	    _limitForScout(colony_size * dimension / 2), 
	    _nb_employed(colony_size/2), 
	    _nb_onlookers(colony_size/2), 
	    _probabilities(0), 
	    _foodSources(0), 
	    _bestSource(), 
	    _lbound(lbound), 
	    _ubound(ubound), 
	    _stop_criteria(stop_criteria), 
	  _cost_function(cost_function)
	{
	  // Initialize our populations
	  _foodSources = new FoodSourceType[_nb_employed];
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _foodSources[i] = FoodSourceType(_dimension);

	  // And the probabilities of their solutions
	  _probabilities = new double[_nb_employed];

	}

	virtual ~Base(void)
	{
	  delete[] _foodSources;
	  delete[] _probabilities;
	}

	void init(void)
	{
	  // Initialize the positions and fitnesses
	  // of the employed bees
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _foodSources[i].init(_lbound, _ubound, _cost_function);

	  _bestSource = _foodSources[0];

	  _epoch = 0;

	  // Keep track of the best solution
	  findBestSource();
	}

	int getEpoch(void)
	{
	  return _epoch;
	}

	void step(void)
	{
	  // Employed bees phase
	  employedPhase();
	  
	  // Onlooker bees phase
	  onlookerPhase();

	  // Scout bees phase
	  scoutPhase();

	  // Memorize the best solution
	  findBestSource();

	  _epoch ++;
	}

	FoodSourceType& getBest()
	{
	  return _bestSource;
	}

	void run(void)
	{
	  while(!_stop_criteria(_bestSource.getFValue(),_epoch))
	    {
	      step();
	    }
	}

	void print(void)
	{
	  std::cout << "Artificial Bee Colony " << std::endl;
	  std::cout << "Food Sources : " << std::endl;
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    std::cout << _foodSources[i] << std::endl;
	  std::cout << "Best source : " << _bestSource << std::endl;

	}
      }; // class Base
    } // namespace algorithm
  } // namespace ABC


  namespace algorithm
  {

    /**
     * ABC algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::ABC::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION> 
    abc(size_t colony_size, size_t dimension,
	const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	const STOP_CRITERIA& stop, const COST_FUNCTION& func)  {
      return popot::ABC::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION>(colony_size, dimension, lbound, ubound, stop, func);
    };

    // The parameters for updating the velocity of the particles
    class SPSO2006_Params
    {
    public:
      static double w(void) { return 1.0/(2.0*log(2.0));};
      static double c(void) { return 0.5 + log(2.0);};
    };
    

    // The particle type, shortcut simplifying the expression of the types in spso2006
    typedef popot::PSO::particle::Particle<> ParticleSPSO;

    /**
     * Builds the SPSO2006 algorithm
     */

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      spso2006(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	       const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function, bool reevaluate_best=false) {

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2006<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      auto init_position_function = popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;
      auto init_velocity_function = popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;
      
      // The confinement method
      auto confine = popot::confinement::confine<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;

      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO, 3, true>;
      //auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO, 3, true>;


      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, confine, init_position_function, init_velocity_function,
					      p, popot::PSO::algorithm::ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION, reevaluate_best); 
      return algo;
    }


    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    spso2006(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));
      return spso2006(swarm_size, dimension, lbound, ubound, stop, cost_function, false);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      stochastic_spso2006(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      return spso2006(swarm_size, dimension, lbound, ubound, stop, cost_function, true);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    stochastic_spso2006(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      return spso2006(swarm_size, dimension, lbound, ubound, stop, cost_function, true);
    }


    /**
     * Builds the SPSO2007 algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      spso2007(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	       const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function, bool reevaluate_best=false) {

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2007<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      auto init_position_function = popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;
      auto init_velocity_function = popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;

      // The confinement method
      auto confine = popot::confinement::confine<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO, 3, true>;
      auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO, 3, true>;

      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, confine, init_position_function, init_velocity_function,
					      p, popot::PSO::algorithm::ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION, reevaluate_best); 
      return algo;
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    spso2007(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      return spso2007(swarm_size, dimension, lbound, ubound, stop, cost_function, false);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    stochastic_spso2007(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      return spso2007(swarm_size, dimension, lbound, ubound, stop, cost_function, true);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      stochastic_spso2007(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      return spso2007(swarm_size, dimension, lbound, ubound, stop, cost_function, true);
    }
    

    /**
     * Builds the SPSO2007 algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      spso2011(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	       const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function, bool reevaluate_best=false) {

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2011<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      auto init_position_function = popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;
      auto init_velocity_function = popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;

      // The confinement method
      auto confine = popot::confinement::confine_spso2011<ParticleSPSO::VECTOR_TYPE, LBOUND_FUNC, UBOUND_FUNC>;
    

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO, 3, true>;
      auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO, 3, true>;

      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, confine, init_position_function, init_velocity_function,
					      p, popot::PSO::algorithm::ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION, reevaluate_best); 

      return algo;
    }


    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    spso2011(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      return spso2011(40, dimension, lbound, ubound, stop, cost_function, false);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
      stochastic_spso2011(size_t swarm_size, size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      return spso2011(swarm_size, dimension, lbound, ubound, stop, cost_function, true);
    }

    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      popot::PSO::algorithm::Base<LBOUND_FUNC, 
				  UBOUND_FUNC, 
				  STOP_CRITERIA,  
				  COST_FUNCTION,
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				  void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
				ParticleSPSO>
    stochastic_spso2011(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION &cost_function) {

      return spso2011(40, dimension, lbound, ubound, stop, cost_function, true);
    }

    /* /\** */
    /*  * Builds the barebone algorithm */
    /*  *\/ */
    /* template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION> */
    /* popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION,  */
    /* 				void(*)(std::vector<ParticleSPSO >&,  */
    /* 					std::vector< typename ParticleSPSO::NeighborhoodType *> &,  */
    /* 					std::map< size_t, std::vector<size_t> > &), */
    /* 				void(*)(ParticleSPSO&), */
    /* 				void(*)(ParticleSPSO&),  */
    /* 				void(*)(ParticleSPSO&), */
    /* 				ParticleSPSO> */
    /* barebone(size_t dimension, */
    /* 			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, */
    /* 			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) { */
    /*   size_t swarm_size = 40; */

    /*   // Particle type */
    /*   ParticleSPSO p; */

    /*   // Position and velocity updates */
    /*   auto position_update = popot::PSO::particle::updatePosition_barebone<ParticleSPSO>; */
    /*   auto velocity_update = popot::PSO::particle::updateVelocity_barebone<ParticleSPSO>; */
     
    /*   // Initialization functions */
    /*   /\* */
    /* 	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void {  */
    /* 	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound); */
    /* 	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound); */
    /* 	}; */
    /*   *\/ */

    /*   // The rule to update the best position */
    /*   auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>; */

    /*   // Topology */
    /*   auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>; */

    /*   auto algo = popot::PSO::algorithm::base(swarm_size, dimension,  */
    /* 					      lbound, ubound, stop, cost_function,  */
    /* 					      topology, position_update, velocity_update, best_position_update,  */
    /* 					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false); */
    /*   return algo; */
    /* } */


    /* /\** */
    /*  * Builds the barebone algorithm */
    /*  *\/ */
    /* template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION> */
    /* popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION,  */
    /* 				void(*)(std::vector<ParticleSPSO >&,  */
    /* 					std::vector< typename ParticleSPSO::NeighborhoodType *> &,  */
    /* 					std::map< size_t, std::vector<size_t> > &), */
    /* 				void(*)(ParticleSPSO&), */
    /* 				void(*)(ParticleSPSO&),  */
    /* 				void(*)(ParticleSPSO&), */
    /* 				ParticleSPSO> */
    /* modified_barebone(size_t dimension, */
    /* 			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, */
    /* 			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) { */
    /*   size_t swarm_size = 40; */

    /*   // Particle type */
    /*   ParticleSPSO p; */

    /*   // Position and velocity updates */
    /*   auto position_update = popot::PSO::particle::updatePosition_barebone<ParticleSPSO>; */
    /*   auto velocity_update = popot::PSO::particle::updateVelocity_modifiedBarebone<ParticleSPSO>; */
     
    /*   // Initialization functions */
    /*   /\* */
    /* 	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void {  */
    /* 	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound); */
    /* 	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound); */
    /* 	}; */
    /*   *\/ */

    /*   // The rule to update the best position */
    /*   auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>; */

    /*   // Topology */
    /*   auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>; */

    /*   auto algo = popot::PSO::algorithm::base(swarm_size, dimension,  */
    /* 					      lbound, ubound, stop, cost_function,  */
    /* 					      topology, position_update, velocity_update, best_position_update,  */
    /* 					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false); */
    /*   return algo; */
    /* } */

  } // namespace algorithm

} // namespace popot

#endif // POPOT_ALGORITHM_H
