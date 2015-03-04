/*   This file is part of popot
 *
 *   Copyright (C) 2011-2012
 *
 *   Author : Jeremy Fix
 *
 *   Contributor :
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *   Contact : Jeremy.Fix@gmail.com
 *
 */

#ifndef POPOT_H
#define POPOT_H

#define VERBOSE_BENCH false

#include "algorithm.h"
#include "benchmark.h"
#include "exceptions.h"
#include "individuals.h"
#include "initializers.h"
#include "maths.h"
#include "mersenne.h"
#include "neighborhood.h"
#include "problems.h"
#include "rng_generators.h"
#include "standard_pso.h"
#include "topology.h"
#include "tools.h"



/*! \mainpage POPulation based Optimization Toolbox
 *
 *  NOT YET UPDATED !!!!!! The 2.0 version introduces major modifications in the way to use the library
 * 
 * \section intro_sec Introduction
 *
 * In this toolbox, you will find various black box optimization algorithms, namely :
 * - Particle Swarm Optimization (with some of its variants), see Engelbrecht
 * - ... more to come 
 *
 * \subsection pso_sec Particle Swarm Optimization
 * Generally speaking, a particle swarm algorithm works by spreading a set of particles
 * in a search space and by moving these particles according to certain rules. The traditional update rules, as defined by
 * Kennedy, involve linear update of the positions and velocities. The update of the velocity
 * depends on the velocity itself (with an inertia factor), on the difference between the best position
 * the particle ever had and its current position (this is called the cognitive component), and the difference between
 * the best position that particles in the neighborhood of the updated particle ever had and its current position
 * (this is called the social component). Different variations have been introduced but the basic concept remains the same.<br/>
 *
 * In the implementation of this library, a swarm is a parametrized object. It can be parametrized by numerical parameters (such as cognitive or social factors)
 * the type of particle, the type of neighborhood, etc... and the problem to solve. The algorithms are defined here as <STRONG>minimization algorithms</STRONG>. Therefore, the problem (see popot::problems namespace for example) provides
 * a function to minimize.
 *
 * Generally, to use the library you need to define :
 * - A problem (see e.g. popot::problems) which provides a function to minimize
 * - A particle type (see e.g. popot::PSO::particle) which takes a problem as a parameter and some parameters (for example the cognitive and social factors)
 * - A topology (see e.g. popot::PSO::topology) which defines the network through which the particles communicate
 * - An algorithm (see e.g. popot::PSO::algorithm) which defines the different steps to update the state of the particles. It requires a particle type, a toplogy and a stopping criteria
 *
 * However, there are already some standard versions of PSO implemented : SPSO2011 and SPSO2006
 * \section sec_problems Problems
 *
 * Different benchmark problems are already coded in the popot::problems namespace :
 * - N-dimensional popot::problems::Ackley function
 * - N-dimensional popot::problems::Griewank function
 * - N-dimensional popot::problems::Sphere function
 * - N-dimensional popot::problems::QuarticNoise function
 * - N-dimensional popot::problems::Rastrigin function
 * - N-dimensional popot::problems::Rosenbrock function
 * - N-dimensional popot::problems::Schwefel1_2 function
 * - N-dimensional popot::problems::Schwefel function
 * - N-dimensional popot::problems::Salomon function
 * - 2 dimensional popot::problems::Dropwave function
 * - ....
 *
 * You can obviously define your own problem. Just check problem.h and the examples to see how one can define its own custom problem.
 *
 * \section sec_particles Particle types
 *
 * Different particles are defined, provided in the popot::PSO::particle namespace :
 * - Standard particle 2006 : popot::PSO::particle::SPSO2006Particle
 * - Standard particle 2011 : popot::PSO::particle::SPSO2011Particle
 * - Barebone particle : popot::PSO::particle::BareboneParticle
 * - Modified barebone particle : popot::PSO::particle::ModifiedBareboneParticle
 * 
 * \subsection non_stochastic_traditional_particle Non stochastic traditional particle
 * This standard PSO relies on the following particle's position and velocity update :<BR>
 * \f{eqnarray*}{
 * v_{k+1} &=& w v_k + c_1 * R_1^T (p^l - p_k) + c_2 * R_2^T (p^g - p_k) \\
 * p_{k+1} &=& p_k + v_{k+1} \\
 * &\mbox{with}& v_k, p_k \in R^N
 * \f}
 * where \f$R_1, R_2\f$ are N-dimensional vectors with uniformely distributed random values in [0;1]
 *
 * The particles can be updated synchronously or asynchronously.
 *
 * \subsection traditional_particle Traditional particle
 *
 * This particle uses the same update rules as the previous particle type except that we reevaluate the fitness
 * of the best position before update the best position. This is usefull in case the fitness is stochastic.
 *
 * \subsection barebone_particle Barebone particle
 *
 * later....
 *
 * \subsection modified_barebone_particle Modified barebone particle
 *
 * later ...
 *
 * \section sec_topologies Topologies
 *
 * The particles communicate with the particles within their neighborhood. The neighborhood is used to define who informs whom, 
 * involved in updating the velocity of the particles. Different topologies are implemented :
 * - Full connectivity (popot::PSO::topology::Full) : the best swarm's position, for each particle, is computed within the \f$N\f$ other particles
 * - Ring connectivity (popot::PSO::topology::Ring) : the best swarm's position, for each particle, is computed within the 2 closest neighbours (left and right) and the best position of the particle itself
 * - VonNeuman (popot::PSO::topology::VonNeuman): the best swarm's position, for each particle, is computed within the 4 closest neighboors (north, east, west, south) and itself. The particles are somehow placed on a 2D grid.
 * - RandomInformants (popot::PSO::topology::RandomInformants) each particle informs K other particles
 * - Probabilistic informants (popot::PSO::topology::AdaptiveRandom) : each particle informs each other particle with a given probability 
 *
 * Each particle hosts a popot::PSO::neighborhood::Neighborhood which is set from the defined topology. A neighborhood can be recomputed for e.g. when there is no improvement of the best position ever found by the algorithm.
 *
 * Shown below are illustrations of three topologies for a swarm with 9 particles (see example-002.cc).
 *
 * \image html "graph_full.png" "Full connectivity" \image html "graph_ring.png" "Ring connectivity" \image html "graph_vonneuman.png" "VonNeuman connectivity"
 *
 * \section sec_algorithms Algorithms
 *
 * So far, a single algorithm is implemented, see popot::PSO::algorithm::Base .
 *
 * \section papers_sec Some relevant papers
 * - Shi, Y. H., Eberhart, R. C., (1998).<STRONG>Parameter Selection in Particle Swarm Optimization</STRONG> , The 7th Annual Conference on Evolutionary Programming, San Diego, USA. (introduced the inertia factor)
 * - Kennedy, J. and Mendes, R. (2002). <STRONG>Population Structure and Particle Swarm Performance</STRONG>
 * - Clerc, M. , Kennedy J. (2002) <STRONG> The Particle Swarm—Explosion, Stability, and Convergence in a Multidimensional Complex Space</STRONG>.
 * - M.E.H. Pedersen (2010), <STRONG>Good Parameters for Particle Swarm Optimization</STRONG>,<EM>Technical Report no. HL1001</EM>
 * - A. Engelbrecht (2010), <STRONG>Heterogeneous Particle Swarm Optimization</STRONG>, <EM>ANTS 2010, LNCS 6234</EM>, pp 191-202
 * - A. Engelbrecht (2010), <STRONG>Fundamentals of computational swarm intelligence</STRONG>, Wiley.
 *
  */

/**
 * \example example-000.cc
 * \example example-001.cc
 * \example example-002.cc
 * \example example-003.cc
 */


#endif // POPOT_H
