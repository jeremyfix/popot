#ifndef POPOT_STANDARD_PSO
#define POPOT_STANDARD_PSO

#include "topology.h"
#include "popot.h"


namespace popot
{
  namespace PSO
  {
    /**
     * SPSO-2006 , Code based on "Standard PSO Descriptions", M. Clerc
     */
    namespace SPSO2006
    {
/*       class Particle_SPSO2006_Params */
/*       { */
/*       public: */
/* 	static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter */
/* 	static double c() { return 0.5 + log(2.0);}   // Best particle position weight */
/*       }; */

/*       template <typename PROBLEM> */
/* 	using Particle = popot::PSO::particle::SPSO2006Particle< PROBLEM, Particle_SPSO2006_Params>; */

/*       /\* template<typename PROBLEM> *\/ */
/*       /\* 	struct Particle *\/ */
/*       /\* 	{ *\/ */
/*       /\* 	  typedef popot::PSO::particle::SPSO2006Particle< PROBLEM, Particle_SPSO2006_Params> Type; *\/ */
/*       /\* 	}; *\/ */

/*       /\** */
/*        * The swarm size is set to 10 + [2 sqrt(S)] */
/*        *\/ */
/*       template<typename PROBLEM> */
/* 	class SwarmSize */
/* 	{ */
/* 	public: */
/* 	  static const int swarm_size; */
/* 	}; */
/*       template<typename PROBLEM> const int SwarmSize<PROBLEM>::swarm_size = 10 + int(2.0 * sqrt(PROBLEM::dimension)); */

/*       /\** */
/*        * Topology of SPSO2006 : Adaptive Random with K=3 */
/*        *\/ */
/*       /\* */
/*       template<typename PROBLEM> */
/* 	struct Topology */
/* 	{ */
/* 	  static const int swarm_size = SwarmSize<PROBLEM>::swarm_size; */
/* 	  typedef popot::PSO::topology::RandomInformants< 3, typename Particle<PROBLEM> > Type; */
/* 	}; */
/*       *\/ */

/*       //TODO : Il faut instancier la topologie quand on connecte les unités, */
/*       // et l'instancier avec la taille  */

/*       template <typename PROBLEM> */
/* 	using Topology = popot::PSO::topology::RandomInformants< 3, Particle<PROBLEM> >; */

/*       class PSO_params  */
/*       { */
/*       public: */
/* 	static bool random_shuffle() { return false;} */
/* 	static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;} */
/*       }; */

/*       template<typename PROBLEM> */
/* 	class Stop_Criteria */
/* 	{ */
/* 	public: */
/* 	  static bool stop(double fitness, int epoch) */
/* 	  { */
/* 	    return  PROBLEM::stop(fitness,epoch); */
/* 	  } */
/* 	}; */
/*       /\* */
/*       template<typename PROBLEM> */
/* 	struct PSO */
/* 	{ */
/* 	  typedef popot::PSO::algorithm::Base<PSO_params, typename Particle<PROBLEM>::Type,typename Topology<PROBLEM>::Type, Stop_Criteria<PROBLEM> > Type; */
/* 	}; */
/* *\/ */
/*       template <typename PROBLEM> */
/* 	using PSO = popot::PSO::algorithm::Base<PSO_params, Particle<PROBLEM>, Topology<PROBLEM>, Stop_Criteria<PROBLEM> >; */


    }


/*     namespace StochasticSPSO2006 */
/*     { */
/*       class Particle_SPSO2006_Params */
/*       { */
/*       public: */
/* 	static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter */
/* 	static double c() { return 0.5 + log(2.0);}   // Best particle position weight */
/*       }; */

/*       template<typename PROBLEM> */
/* 	struct Particle */
/* 	{ */
/* 	  typedef popot::PSO::particle::StochasticSPSO2006Particle< PROBLEM, Particle_SPSO2006_Params> Type; */
/* 	}; */

/*       /\** */
/*        * The swarm size is set to 10 + [2 sqrt(S)] */
/*        *\/ */
/*       template<typename PROBLEM> */
/* 	class SwarmSize */
/* 	{ */
/* 	public: */
/* 	  static const int swarm_size; */
/* 	}; */
/*       template<typename PROBLEM> const int SwarmSize<PROBLEM>::swarm_size = 10 + int(2.0 * sqrt(PROBLEM::nb_parameters)); */

/*       /\** */
/*        * Topology of SPSO2006 : Adaptive Random with K=3 */
/*        *\/ */
/*       template<typename PROBLEM> */
/* 	struct Topology */
/* 	{ */
/* 	  static const int swarm_size = SwarmSize<PROBLEM>::swarm_size; */
/* 	  typedef popot::PSO::topology::RandomInformants< swarm_size , 3, typename Particle<PROBLEM>::Type> Type; */
/* 	}; */

/*       class PSO_params */
/*       { */
/*       public: */
/* 	static bool random_shuffle() { return false;} */
/* 	static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;} */
/*       }; */
/*       template<typename PROBLEM> */
/*       class Stop_Criteria */
/*       { */
/*       public: */
/* 	static bool stop(double fitness, int epoch) */
/* 	{ */
/* 	  return  PROBLEM::stop(fitness,epoch); */
/* 	} */
/*       }; */

/*       template<typename PROBLEM> */
/* 	struct PSO */
/* 	{ */
/* 	  typedef popot::PSO::algorithm::Base<PSO_params, typename Particle<PROBLEM>::Type,typename Topology<PROBLEM>::Type, Stop_Criteria<PROBLEM> > Type; */
/* 	}; */

/*     } */


/*     /\** */
/*      * SPSO-2007 , Code based on "Standard PSO Descriptions", M. Clerc */
/*      *\/ */
/*     namespace SPSO2007 */
/*     { */

/*       class Particle_SPSO2007_Params */
/*       { */
/*       public: */
/* 	static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter */
/* 	static double c() { return 0.5 + log(2.0);}   // Best particle position weight */
/*       }; */
/*       template<typename PROBLEM> */
/* 	struct Particle */
/* 	{ */
/* 	  typedef popot::PSO::particle::SPSO2007Particle< PROBLEM, Particle_SPSO2007_Params> Type; */
/* 	}; */


/*       /\** */
/*        * The swarm size is set to 10 + [2 sqrt(S)] */
/*        *\/ */
/*       template<typename PROBLEM> */
/* 	class SwarmSize */
/* 	{ */
/* 	public: */
/* 	  static const int swarm_size; */
/* 	}; */
/*       template<typename PROBLEM> const int SwarmSize<PROBLEM>::swarm_size = 10 + int(2.0 * sqrt(PROBLEM::nb_parameters)); */

/*       /\** */
/*        * Topology of SPSO2006 : Adaptive Random with K=3 */
/*        *\/ */
/*       template<typename PROBLEM> */
/* 	struct Topology */
/* 	{ */
/* 	  static const int swarm_size = SwarmSize<PROBLEM>::swarm_size; */
/* 	  typedef popot::PSO::topology::AdaptiveRandom< swarm_size , 3, typename Particle<PROBLEM>::Type> Type; */
/* 	}; */


/*       class PSO_params */
/*       { */
/*       public: */
/* 	static bool random_shuffle() { return false;} */
/* 	static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;} */
/*       }; */
/*       template<typename PROBLEM> */
/*       class Stop_Criteria */
/*       { */
/*       public: */
/* 	static bool stop(double fitness, int epoch) */
/* 	{ */
/* 	  return  PROBLEM::stop(fitness,epoch); */
/* 	} */
/*       }; */

/*       template<typename PROBLEM> */
/* 	struct PSO */
/* 	{ */
/* 	  typedef popot::PSO::algorithm::Base<PSO_params, typename Particle<PROBLEM>::Type,typename Topology<PROBLEM>::Type, Stop_Criteria<PROBLEM> > Type; */
/* 	}; */



/*     } */

/*     /\** */
/*      * SPSO-2011 , Code based on "Standard PSO Descriptions", M. Clerc */
/*      *\/ */
/*     namespace SPSO2011 */
/*     { */
/*       class Particle_SPSO2011_Params */
/*       { */
/*       public: */
/* 	static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter */
/* 	static double c() { return 0.5 + log(2.0);}   // Best particle position weight */
/*       }; */
/*       template<typename PROBLEM> */
/* 	struct Particle */
/* 	{ */
/* 	  typedef popot::PSO::particle::SPSO2011Particle< PROBLEM, Particle_SPSO2011_Params> Type; */
/* 	}; */

/*       template<typename PROBLEM> */
/* 	struct Topology */
/* 	{ */
/* 	  static const int swarm_size = 40; */
/* 	  typedef popot::PSO::topology::AdaptiveRandom<swarm_size , 3, typename Particle<PROBLEM>::Type> Type; */
/* 	}; */


/*       class PSO_params */
/*       { */
/*       public: */
/* 	static bool random_shuffle() { return false;} */
/* 	static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;} */
/*       }; */
/*       template<typename PROBLEM> */
/*       class Stop_Criteria */
/*       { */
/*       public: */
/* 	static bool stop(double fitness, int epoch) */
/* 	{ */
/* 	  return  PROBLEM::stop(fitness,epoch); */
/* 	} */
/*       }; */

/*       template<typename PROBLEM> */
/* 	struct PSO */
/* 	{ */
/* 	  typedef popot::PSO::algorithm::Base<PSO_params, typename Particle<PROBLEM>::Type,typename Topology<PROBLEM>::Type, Stop_Criteria<PROBLEM> > Type; */
/* 	}; */
/*     } */


/*     /\** */
/*      * Stochastic SPSO-2011 , Code based on "Standard PSO Descriptions", M. Clerc */
/*      *\/ */
/*     namespace StochasticSPSO2011 */
/*     { */
/*       class Particle_SPSO2011_Params */
/*       { */
/*       public: */
/* 	static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter */
/* 	static double c() { return 0.5 + log(2.0);}   // Best particle position weight */
/*       }; */
/*       template<typename PROBLEM> */
/* 	struct Particle */
/* 	{ */
/* 	  typedef popot::PSO::particle::StochasticSPSO2011Particle< PROBLEM, Particle_SPSO2011_Params> Type; */
/* 	}; */

/*       template<typename PROBLEM> */
/* 	struct Topology */
/* 	{ */
/* 	  static const int swarm_size = 40; */
/* 	  //typedef popot::PSO::topology::FixedRandomInformants< swarm_size , 5, typename Particle<PROBLEM>::Type> Type; */
/* 	  typedef popot::PSO::topology::AdaptiveRandom<swarm_size , 3, typename Particle<PROBLEM>::Type> Type; */
/* 	}; */


/*       class PSO_params */
/*       { */
/*       public: */
/* 	static bool random_shuffle() { return false;} */
/* 	static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;} */
/*       }; */
/*       template<typename PROBLEM> */
/*       class Stop_Criteria */
/*       { */
/*       public: */
/* 	static bool stop(double fitness, int epoch) */
/* 	{ */
/* 	  return  PROBLEM::stop(fitness,epoch); */
/* 	} */
/*       }; */

/*       template<typename PROBLEM> */
/* 	struct PSO */
/* 	{ */
/* 	  typedef popot::PSO::algorithm::Base<PSO_params, typename Particle<PROBLEM>::Type,typename Topology<PROBLEM>::Type, Stop_Criteria<PROBLEM> > Type; */
/* 	}; */
/*     } */



  }
}



#endif
