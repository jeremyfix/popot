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

#ifndef POPOT_PROBLEMS_RL_H
#define POPOT_PROBLEMS_RL_H

#include <tuple>
#include <cmath>
#include <fstream>

namespace popot {

  namespace problems {
    /**
     * Problems in reinforcement learning
     **/
    namespace rl {

      /**
       * The definition of the mountain car problem
       **/
      namespace mountain_car {

	struct Params {
	  // The boundaries of the problem
	  static std::pair<double, double> x_range() { return std::make_pair(-1.2, 0.5);};
	  static std::pair<double, double> dx_range() { return std::make_pair(-0.07, 0.07);};

	  // The ranges for initializing the car
	  static std::pair<double, double> x0_range() { return std::make_pair(-0.75, -0.25);};
	  static std::pair<double, double> dx0_range() { return std::make_pair(-0.02, 0.02);};

	  // Maximal number of iterations allowed for a trial
	  static unsigned int max_length_episode() { return 1500;};
	};

	template<typename PARAMS>
	struct State {
	  double _x, _dx;

	  State() {
	    _x = popot::math::uniform_random(PARAMS::x0_range().first, PARAMS::x0_range().second);
	    _dx = popot::math::uniform_random(PARAMS::dx0_range().first, PARAMS::dx0_range().second);
	  }
	};

	struct Action {
	  double _f;
	};

	template<typename TPOLICY, typename PARAMS>
	struct Simulator {
	  static void transition(State<PARAMS>& s, Action& a, double & rew, int &eoe) {
	    double aa = a._f;

	    double x_t = s._x;
	    double dx_t = s._dx;

	    // Update the velocity
	    double dx_t_1 = dx_t + 1e-3 * (aa - 2.5 * cos(3.0 * x_t));

	    // Bound it
	    if(dx_t_1 <= PARAMS::dx_range().first)
	      dx_t_1 = PARAMS::dx_range().first;
	    else if(dx_t_1 >= PARAMS::dx_range().second)
	      dx_t_1 = PARAMS::dx_range().second;

	    // Update the position
	    double x_t_1 = x_t + dx_t_1;

	    // Bound it
	    if(x_t_1 <= PARAMS::x_range().first) {
	      x_t_1 = PARAMS::x_range().first;
	      dx_t_1 = 0.0;
	    }

	    if(x_t_1 >= PARAMS::x_range().second) {
	      x_t_1 = PARAMS::x_range().second;
	      dx_t_1 = 0.0;
	    }

	    // Set the new position and velocity
	    s._x = x_t_1;
	    s._dx = dx_t_1;

	    if(x_t_1 >= PARAMS::x_range().second) {
	      rew = 0;
	      eoe = 1;
	    }
	    else {
	      rew = -1.0;
	      eoe = 0;
	    }
	  }


	  static double evaluate(double * params) {
	    return evaluaten(params, 1);
	  }

	  static double evaluaten(double * params, int nb_trajectories) {
	    int epoch;
	    int eoe;
	    double rew = 0.0;
	    double rew_t;
	    double rew_cumul = 0.0;
	    for(int i = 0 ; i < nb_trajectories; ++i) {
	      State<PARAMS> s;
	      Action a;
	      epoch = 0;
	      eoe = 0;
	      rew_cumul = 0;
	      while((epoch < PARAMS::max_length_episode()) && !eoe) {
		TPOLICY::sample_action(s, params, a);
		transition(s, a, rew_t, eoe);
		rew_cumul += rew_t;
		++epoch;
	      }
	      rew += rew_cumul;
	    }
	    return rew / double(nb_trajectories);
	  }

	  // This functions test a policy on one trial and prints the state, etc.. in filename
	  static void print_policy(std::string filename, double * params)
	  {
	    std::ofstream outfile;

	    State<PARAMS> s;
	    Action a;

	    outfile.open(filename.c_str());
	    int epoch = 0;
	    int eoe=0;
	    double rew_t;

	    while( (epoch < PARAMS::max_length_episode()) && !eoe) {
	      TPOLICY::sample_action(s, params, a);
	      transition(s, a, rew_t, eoe);
	      ++epoch;
	      outfile << epoch << "\t"
		      << s._x << "\t"
		      << s._dx << "\t"
		      << a._f << "\t"
		      << rew_t << std::endl;
	    }
	    outfile.close();

	    std::cout << "Saved a policy in " << filename << std::endl;
	    std::cout << "Columns must be understood as : epoch x dx action reward" << std::endl;
	  }

	  // This function evaluates the average number steps to reach the goal
	  // averaged over NB_POLICY_EVALUATIONS
	  static double test_policy(double * params, int nb_trajectories) {
	    double total_steps_to_goal = 0.0;
	    int epoch, eoe;
	    double rew_t;
	    for(int i = 0 ; i < nb_trajectories ; ++i)
	      {
		epoch = 0;
		eoe=0;
		State<PARAMS> s;
		Action a;
		while( (epoch < PARAMS::max_length_episode()) && !eoe)
		  {
		    TPOLICY::sample_action(s, params, a);
		    transition(s, a, rew_t, eoe);
		    ++epoch;
		  }
		total_steps_to_goal += epoch;
	      }
	    return total_steps_to_goal/double(nb_trajectories);
	  }	  
	  
	  // This function evaluates the average number of steps to reach the goal
	  // from uniformely distributed starting positions (x, dx)
	  static void test_length_all_positions(std::string filename, double * params, int nb_trajectories) {
	    State<PARAMS> s;
	    Action a;

	    double rew_t;
	    int eoe;
	    int epoch = 0;
	    int N = 50;

	    double init_x, init_dx;
	    double total_length_episodes = 0;
	    std::ofstream outfile(filename.c_str());

	    double x_min, x_max, dx_min, dx_max;
	    std::tie(x_min, x_max) = PARAMS::x_range();
	    std::tie(dx_min, dx_max) = PARAMS::dx_range();
	    
	    for(int i = 0 ; i < N ; ++i)
	      {
		init_x = x_min + (x_max - x_min) *double(i)/double(N-1);
		for(int j = 0 ; j < N ; ++j)
		  {
		    init_dx = dx_min + (dx_max - dx_min) *double(j)/double(N-1);
		    total_length_episodes = 0;
		    for(int k = 0 ; k< nb_trajectories ; ++k)
		      {
			s._x = init_x;
			s._dx = init_dx;
			epoch = 0;
			eoe = 0;

			while((epoch < PARAMS::max_length_episode()) && !eoe)
			  {
			    TPOLICY::sample_action(s, params, a);
			    transition(s, a, rew_t, eoe);
			    ++epoch;
			  }
			total_length_episodes += epoch;
		      }
		    outfile << init_x << " " << init_dx << " " << total_length_episodes/double(nb_trajectories) << std::endl;
		  }
	      }
	    outfile.close();
	  }


	};


	struct RBFParams {
	  static int nb_centers() { return 3;}
	  static double sigma() { return 0.3;}
	};

	/**
	 * Policy : RBF with 9 gaussians + a constant term  per action 
	 * Selection is done with a Gibbs schema over         
	 * these approximations of Q(s,a)             
	 **/
	template<typename SPARAMS, typename PPARAMS>
	struct RBFPolicy {

	  static int dimension() {
	    return 3 * (PPARAMS::nb_centers() * PPARAMS::nb_centers() + 1);
	  }

	  static double lbound(size_t index) { return -100;}
	  static double ubound(size_t index) { return  100;}
	  


	  // Compute the probabilities to select each of the 3 actions
	  // These probabilites, returned in the proba array are normalized, they sum to 1.0
	  static double compute_probability(State<SPARAMS> &s, double * params, double * proba)
	  {
	    double mu_x, mu_dx;
	    double phi_s_a;
	    double scaled_x, scaled_dx;
	    double sum_proba = 0.0;
	    int index_param = 0;
	    double temperature = 1.0;

	    double x_min, x_max, dx_min, dx_max;
	    std::tie(x_min, x_max) = SPARAMS::x_range();
	    std::tie(dx_min, dx_max) = SPARAMS::dx_range();

	    scaled_x = 1.0 - (x_max - s._x)/(x_max - x_min);
	    scaled_dx = 1.0 - (dx_max - s._dx)/(dx_max- dx_min);

	    for(int i = 0 ; i < 3 ; ++i)  {
	      // The constant term
	      proba[i] = params[index_param];
	      ++index_param;
	      for(int j_x = 0 ; j_x < PPARAMS::nb_centers() ; ++j_x) {
		// The gaussians on x are centered on [0, 0.5 , 1.0]
		mu_x =  double(j_x)/ double(PPARAMS::nb_centers() - 1.0);
		for(int j_dx = 0 ; j_dx < PPARAMS::nb_centers() ; ++j_dx) {
		  // The gaussians on dx are centered on [0, 0.5 , 1.0]
		  mu_dx = double(j_dx) / double(PPARAMS::nb_centers() - 1.0);
		  phi_s_a = exp(-((mu_x - scaled_x)*(mu_x - scaled_x) + (mu_dx - scaled_dx)*(mu_dx - scaled_dx))/(2.0*PPARAMS::sigma()*PPARAMS::sigma()));
		  proba[i] = proba[i] + params[index_param]*phi_s_a ;
		  index_param++;
		}
	      }
	      proba[i] = exp(proba[i]/temperature);
	      sum_proba = sum_proba + proba[i];
	    }
	    
	    return sum_proba;
	  }
	  
	  static void sample_action(State<SPARAMS> &s, double * params, Action &a)
	  {
	    double proba[3];

	    // Compute the probabilities to select each action
	    // normalized in 0, 1.0
	    double sum_proba = compute_probability(s, params, proba);

	    double rand_val_action = sum_proba * popot::math::uniform_random(0.0, 1.0);
	    double str;
	    if(rand_val_action <= proba[0])
	      str = -1.0;
	    else if(rand_val_action <= proba[0]+proba[1])
	      str = 0.0;
	    else
	      str = 1.0;
	    a._f = str;

	  }
	};




      }

      /**
       * The definition of the inverted pendulum (cart-pole) problem
       **/
      namespace inverted_pendulum {
	struct Params {
	  static std::pair<double, double> theta_range() { return std::make_pair(-1e-1, 1e-1);}
	  static std::pair<double, double> dtheta_range() { return std::make_pair(-1e-1, 1e-1);}
	  static unsigned int max_length_episode() { return 3000;}
	};

	struct Pendulum {
	  static double g(void)        {return 9.8;}
	  static double m(void)        {return 2.0;}
	  static double M(void)        {return 8.0;}
	  static double l(void)        {return 0.5;}
	  static double a(void)        {return 1.0/(m()+M());}
	  static double strength(void) {return 50.0;}
	  static double tau(void)      {return  0.1;}
	  static double aml(void)      {return a()*m()*l();}
	  static double actionNoise(void) { return 0.1;} 
	};

	template<typename PARAMS>
	struct State {
	  double _theta, _dtheta;
	  
	  State() {
	    _theta = popot::math::uniform_random(PARAMS::theta_range().first, PARAMS::theta_range().second);
	    _dtheta = popot::math::uniform_random(PARAMS::dtheta_range().first, PARAMS::dtheta_range().second);
	  }
	};
       
	struct Action {
	  double _f;
	};

	template<typename TPOLICY, typename StateParams, typename PendulumParams>
	class Simulator
	{
	public:
	  static void transition(State<StateParams> &s, Action& a, double &rew, int &eoe) {
	    double aa = a._f;
	    aa = aa + PendulumParams::actionNoise()*popot::math::uniform_random(-1., 1.);
	    aa = PendulumParams::strength() * aa;

	    double theta = s._theta;
	    double dtheta = s._dtheta;

	    double cphi = cos(theta);
	    double acc = ( PendulumParams::g()*sin(theta)
			   - 0.5*PendulumParams::aml()*sin(2*theta)*dtheta*dtheta
			   - PendulumParams::a()*cphi*aa )
	      / ( 4*PendulumParams::l()/3.0 - PendulumParams::aml()*cphi*cphi );

	    s._theta = theta + dtheta*PendulumParams::tau();
	    s._dtheta = dtheta + acc * PendulumParams::tau();
	    
	    if(fabs(s._theta)>= M_PI/2.0) {
	      rew = -1.0;
	      eoe = 1;
	    }
	    else
	      {
		rew = 0.0;
		eoe = 0;
	      }
	  }

	  static double evaluate(double * params) {
	    return evaluaten(params, 1);
	  }

	  static double evaluaten(double * params, int nb_trajectories) {
	    int epoch;
	    int eoe;
	    double rew=0.0;
	    double rew_t;
	    double rew_cumul=0;

	    for(int i = 0 ; i < nb_trajectories ; ++i)
	      {
		epoch = 0;
		eoe = 0;
		rew_cumul = 0;

		State<StateParams> s;
		Action a;

		while( (epoch < StateParams::max_length_episode()) && !eoe)
		  {
		    TPOLICY::sample_action(s, params, a);
		    transition(s, a, rew_t, eoe);
		    rew_cumul += rew_t;
		    epoch++;
		  }
		rew += rew_cumul/double(epoch);
	      }

	    return rew/double(nb_trajectories);
	  }

	  static void print_policy(std::string filename, double * params)
	  {
	    std::ofstream outfile;

	    State<StateParams> s;
	    Action a;

	    outfile.open(filename.c_str());
	    int epoch = 0;
	    int eoe=0;
	    double rew_t;

	    while( (epoch < StateParams::max_length_episode()) && !eoe)
	      {
		// Sample an action according to the policy
		TPOLICY::sample_action(s, params, a);
		// Make the transition
	        transition(s, a, rew_t, eoe);
		epoch++;
		outfile << epoch << "\t"
			<< s._theta << "\t"
			<< s._dtheta << "\t"
			<< a._f << "\t"
			<< rew_t << std::endl;
	      }
	    outfile.close();

	    std::cout << "I saved a policy in " << filename << "; columns must be read as : epoch theta dtheta action reward" << std::endl;
	  }



	  // This function evaluates the average number of balancing steps of a policy
	  // averaged over nb_trajectories
	  static double test_policy(double * params, int nb_trajectories)
	  {
	    double total_steps_balancing = 0.0;
	    int epoch, eoe;
	    double rew_t;

	    for(int i = 0 ; i < nb_trajectories ; ++i)
	      {
		epoch = 0;
		eoe=0;
		// Initialize randomly the pendulum
		State<StateParams> s;
		Action a;

		while( (epoch < StateParams::max_length_episode()) && !eoe)
		  {
		    // Sample an action according to the policy
		    TPOLICY::sample_action(s, params, a);
		    // Make the transition
		    transition(s, a, rew_t, eoe);
		    epoch++;
		  }
		total_steps_balancing = total_steps_balancing + epoch;
	      }
	    return total_steps_balancing/double(nb_trajectories);
	  }
	};

	/**
	 * Policy whose actions is drawn from a normal distribution
	 * of mean t0 * theta + t1 * dtheta 
	 * and variance 1.0 / (1.0 + exp(-t2))
	 **/
	template<typename STATE>
	class NormalPolicy
	{
	public:
	  static int dimension() {
	    return 3;
	  }
	  static double lbound(size_t index) { return -100;}
	  static double ubound(size_t index) { return  100;}
	  

	  static void sample_action(STATE &s, double * params, Action &a) {
	    double sigma = 1.0 / (1.0 + exp(-params[2]));
	    double mean = params[0] * s._theta + params[1] * s._dtheta;
	    double str = popot::math::normal(mean,sqrt(sigma));

	    a._f = str;
	  }

	  static void sample_random_action(STATE & s, Action &a) {
	    a._f = popot::math::uniform_random(-1.0, 1.0);
	  }
	};


	/**
	 *  Policy : RBF with 9 gaussians + a constant term  per action
	 *  Selection is done with a Gibbs schema over  
         *  these approximations of Q(s,a)     
	 **/

	struct RBFParams{
	  static int nb_centers() { return 3;}
	  static double sigma() { return 1.0;}
	};

	template<typename STATE, typename PARAMS>
        struct RBFPolicy
	{
	  static int dimension() {
	    return 3 * (PARAMS::nb_centers() * PARAMS::nb_centers() + 1);
	  }
	  static double lbound(size_t index) { return -100;}
	  static double ubound(size_t index) { return  100;}

	  // Compute the probabilities to select each of the 3 actions
	  // These probabilites, returned in the proba array are normalized, they sum to 1.0
	  static void compute_probability(STATE &s, double * params, double * proba)
	  {        
	    double mu_theta, mu_dtheta;
	    double phi_s_a;
	    int index_param=0;
	    double sum_proba = 0.0;

	    double theta = s._theta;
	    double dtheta = s._dtheta;

	    double temperature = 1.0;

	    for(int i = 0 ; i < 3 ; ++i)
	      {
		// The constant term
		proba[i] = params[index_param];
		index_param++;
		for(int j_theta = 0 ; j_theta < PARAMS::nb_centers() ; ++j_theta) {
		  // The gaussians on theta are centered on [-Pi/4, 0 , Pi/4]
		  mu_theta = -M_PI_4 + double(j_theta) * M_PI_2 /(PARAMS::nb_centers()-1.);
		  for(int j_dtheta = 0 ; j_dtheta < PARAMS::nb_centers() ; ++j_dtheta)  {
		    // The gaussians on dtheta are centered on [-1 , 0, 1]
		    mu_dtheta = -1.0 + double(j_dtheta) * 2.0 / (PARAMS::nb_centers()-1.);

		    // Compute the basis function
		    phi_s_a = exp(-( (mu_theta - theta)*(mu_theta - theta) + (mu_dtheta - dtheta)*(mu_dtheta - dtheta) )/(2.0*PARAMS::sigma() * PARAMS::sigma()));
		    proba[i] = proba[i] + params[index_param] * phi_s_a;
		    index_param++;
		  }
		}
		proba[i] = exp(proba[i]/temperature);
		sum_proba = sum_proba + proba[i];
	      }

	    // Normalize the probabilities
	    for(int i = 0 ; i < 3 ; ++i)
	      proba[i] = proba[i] / sum_proba;
	  }

	  static void sample_action(STATE &s, double * params, Action &a)
	  {
	    double proba[3];
	    double str,rand_val_action;

	    // Compute the probabilities to select each action
	    // normalized in 0, 1.0
	    compute_probability(s, params, proba);

	    // Toss a random number in [0,1.0]
	    rand_val_action = popot::math::uniform_random(0., 1.0);
	    if(rand_val_action <= proba[0])
	      str = -1.0;
	    else if(rand_val_action <= proba[0]+proba[1])
	      str = 0.0;
	    else
	      str = 1.0;

	    // Set the strength of the action
	    a._f = str;
	  }
	};



      }

      /**
       * Acrobat problem. With the parameters given below, the problem is deterministic
       **/
      namespace acrobat {

	struct Params {
	  // The bounds for theta1 and theta2 must cover a range of 2 pi
	  static std::pair<double, double> theta1_range() { return std::make_pair(-M_PI_2, 3.0 * M_PI_2);}
	  static std::pair<double, double> theta2_range() { return std::make_pair(-M_PI, M_PI);}
	  // The bounds for the speeds
	  static std::pair<double, double> dtheta1_range() { return std::make_pair(-9.0 * M_PI, 9.0 * M_PI);}
	  static std::pair<double, double> dtheta2_range() { return std::make_pair(-4.0 * M_PI, 4.0 * M_PI);}

	  static double delta_t() { return 0.01;}
	  static int max_length_episode() { return 20.0 / delta_t();}
	};

	struct Acrobat {
	  static double strNoise(void)  { return 0.0;}

	  static double gravity(void) { return 9.8;}

	  static double m1(void) { return 1.0;}
	  static double l1(void) { return 1.0;}
	  static double lc1(void) { return l1()/2.0;}
	  static double I1(void) { return 1.0/12.0 * m1() * l1()* l1();}

	  static double m2(void) { return 1.0;}
	  static double l2(void) { return 1.0;}
	  static double lc2(void) { return l2()/2.0;}
	  static double I2(void) { return 1.0/12.0 * m2() * l2()* l2();}

	  static double m1lc1(void) { return m1() * lc1();}
	  static double m2l1(void) { return m2() * l1();}
	  static double m2lc2(void) { return m2() * lc2();}
	  static double m1lc1_2(void) { return m1()*lc1()*lc1();}
	  static double m2l1_2(void) { return m2() * l1() * l1();}
	  static double m2lc2_2(void) { return m2() * lc2() * lc2();}
	  static double m2l1lc2(void) { return m2() * l1() * lc2();}

	  static double max_torque(void) { return 2.0;}
	};

	struct GoalState {
	  static std::pair<double, double> theta1_range() { return std::make_pair(M_PI_2-M_PI/4.0, M_PI_2 + M_PI/4.0);} 
	  static std::pair<double, double> theta2_range() { return std::make_pair(-M_PI/4.0,M_PI/4.0);}

	  static std::pair<double, double> dtheta1_range() { return std::make_pair(-M_PI_2, M_PI_2);}
	  static std::pair<double, double> dtheta2_range() { return std::make_pair(-M_PI_2, M_PI_2);}
	};
	
	static void bound(double& var, std::pair<double, double> bounds, bool angular=false) {
	  if(!angular) {
	    if(var < bounds.first) 
	      var = bounds.first;
	    else if(var > bounds.second)
	      var = bounds.second;
	  }
	  else {
	    if(var < bounds.first)
	      while(var < bounds.first)
		var += 2.0 * M_PI;
	    else if(var > bounds.second)
	      while(var > bounds.second)
		var -= 2.0 * M_PI;
	  }
	}

	static bool is_in_bounds(double var, std::pair<double, double> bounds) {
	  return var >= bounds.first && var <= bounds.second;
	}


	template<typename PARAMS>
	struct State {
	  double _theta1, _theta2, _dtheta1, _dtheta2;

	  State() {
	    // The pendulum is initialized pointing straight down, at rest
	    _theta1 = 3.0 * M_PI_2;
	    _theta2 = 0.0;
	    _dtheta1 = 0.0;
	    _dtheta2 = 0.0;
	  }
	};

	struct Action {
	  double _f;
	};


	template<typename TPOLICY, typename StateParams, typename GoalStateParams, typename AcrobatParams>
	struct Simulator {
	  
	  static int func (double t, const double y[], double dydt[], void * params) {
	    double torque = *((double*)params);
	    double theta_1 = y[0];
	    double dtheta_1 = y[1];
	    double theta_2 = y[2];
	    double dtheta_2 = y[3];

	    double ddtheta_1;
	    double ddtheta_2;

	    double t1 = 0.0;
	    double t2 = torque;

	    double c1 = cos(theta_1);
	    double c2 = cos(theta_2);
	    double c12 = cos(theta_1 + theta_2);
	    double s2 = sin(theta_2);

	    double d11 = AcrobatParams::m1lc1_2() + AcrobatParams::m2l1_2() + AcrobatParams::m2lc2_2() + 2.0 * AcrobatParams::m2l1lc2() * c2 + AcrobatParams::I1() + AcrobatParams::I2();
	    double d12 = AcrobatParams::m2lc2_2() + AcrobatParams::m2l1lc2() * c2 + AcrobatParams::I2();
	    double d21 = AcrobatParams::m2lc2_2() + AcrobatParams::m2l1lc2() * c2 + AcrobatParams::I2();
	    double d22 = AcrobatParams::m2lc2_2() + AcrobatParams::I2();


	    double h1 = -AcrobatParams::m2l1lc2() * s2 * dtheta_2 *dtheta_2 - 2.0 * AcrobatParams::m2l1lc2() * s2 * dtheta_2 * dtheta_1;
	    double h2 = AcrobatParams::m2l1lc2() * s2 * dtheta_1 * dtheta_1;

	    double phi1 = (AcrobatParams::m1lc1() +AcrobatParams::m2l1()) * AcrobatParams::gravity() * c1 + AcrobatParams::m2lc2() * AcrobatParams::gravity() * c12;
	    double phi2 = AcrobatParams::m2lc2() * AcrobatParams::gravity() * c12;


	    double y1 = t1 - h1 - phi1;
	    double y2 = t2 - h2 - phi2;
	    double det = d11 * d22 - d21 * d12;

	    // The evolution equations are defined as :
	    // D * [q1''; q2''] + [h1+phi1 ; h2+phi2] = [t1 ; t2] (t1 = 0; t2 = torque)
	    // <=>  D * [ q1'' ; q2''] = [y1 ; y2]
	    // with :
	    // y1 = t1 - h1 - phi1;
	    // y2 = t2 - h2 - phi2;

	    // Therefore, we get the accelarations by solving the cramer system:
	    // q1'' = 1/det(D) * det[ y1 , d12 ; y2 , d22]
	    //      = 1/(d11 d22 - d21 d12) * ( y1*d22 - y2 d12)
	    // q2'' = 1/det(D) * det[d11 , y1 ; d21 , y2]
	    //      = 1/(d11 d22 - d21 d12) * (d11 * y2 - d21 * y1)

	    if ( fabs(det) < 1e-20 ) { 
	      printf("WARNING : Acrobot singularity\n");
	      return 0;
	    }
	      
	    ddtheta_1 = (1.0/det) * (y1 * d22 - y2 * d12);
	    ddtheta_2 = (1.0/det) * (d11 * y2 - d21 * y1);
	      
	    dydt[0] = dtheta_1; // Dtheta1
	    dydt[1] = ddtheta_1; // DDTheta1
	    dydt[2] = dtheta_2; // Dtheta2
	    dydt[3] = ddtheta_2; // DDTheta2

	    return 1;

	  }

	  static double end_point(State<StateParams> &s)  {
	    return AcrobatParams::l1() * sin(s._theta1) + AcrobatParams::l2() * sin(s._theta1+s._theta2);
	  }
	  

	  static void transition(State<StateParams> &s, Action& a, double &rew, int &eoe) {
	    double torque = a._f + AcrobatParams::strNoise() * popot::math::uniform_random(-1, 1);
	    double y[4] = {s._theta1 , s._dtheta1 , s._theta2, s._dtheta2};
	    double dy[4] = {0.0,0.0,0.0,0.0};
	    func (0.0, y, dy, &torque);
	    
	    s._theta1 += StateParams::delta_t() * dy[0];
	    s._dtheta1 += StateParams::delta_t() * dy[1];
	    s._theta2 += StateParams::delta_t() * dy[2];
	    s._dtheta2 += StateParams::delta_t() * dy[3];
	  
	    bound(s._theta1, StateParams::theta1_range(), true);
	    bound(s._theta2, StateParams::theta2_range(), true);

	    bound(s._dtheta1, StateParams::dtheta1_range(), false);
	    bound(s._dtheta2, StateParams::dtheta2_range(), false);

	    if(is_in_bounds(s._theta1, GoalStateParams::theta1_range()) &&
	       is_in_bounds(s._theta2, GoalStateParams::theta2_range()) &&
	       is_in_bounds(s._dtheta1, GoalStateParams::dtheta1_range()) &&
	       is_in_bounds(s._dtheta2, GoalStateParams::dtheta2_range()))
	      rew = 1.0;
	    else
	      rew = 0.0;
	  }
	  
	  static double evaluate(double * params) {
	    return evaluaten(params, 1);
	  }
	  
	  static double evaluaten(double * params, int nb_trajectories) {
	    int epoch;
	    int eoe;
	    double rew=0.0;
	    double rew_t;
	    double rew_cumul=0;

	    for(int i = 0 ; i < nb_trajectories ; ++i)  {
	      epoch = 0;
	      eoe = 0;
	      rew_cumul = 0;
	      State<StateParams> s;
	      Action a;
	      while( (epoch < StateParams::max_length_episode()) && !eoe) {
		TPOLICY::sample_action(s, params, a);
		transition(s, a, rew_t, eoe);
		rew_cumul += rew_t;
		++epoch;
	      }
	      rew += rew_cumul;
	    }
	    
	    return rew/double(nb_trajectories);
	    
	  }

	  static double test_policy(double* params, int nb_trajectories) {
	    return evaluaten(params, nb_trajectories);
	  }
	  
	  // This functions test a policy on one trial and prints the state, etc.. in filename
	  static void print_policy(std::string filename, double* params) {
	    std::ofstream outfile;

	    State<StateParams> s;
	    Action a;

	    outfile.open(filename.c_str());
	    int epoch = 0;
	    int eoe=0;
	    double rew_t;
	    while( (epoch < StateParams::max_length_episode()) && !eoe)
	      {
		// Sample an action according to the random policy
		TPOLICY::sample_action(s, params, a);
		// Make the transition
		transition(s, a, rew_t, eoe);
		++epoch;
		outfile << epoch << "\t"
			<< s._theta1 << "\t"
			<< s._theta2 << "\t"
			<< s._dtheta1 << "\t"
			<< s._dtheta2 << "\t"
			<< a._f << "\t"
			<< end_point(s) << "\t"
			<< rew_t << std::endl;
	      }
	    outfile.close();

	    std::cout << "Saved a policy in " << filename << std::endl;
	    std::cout << "Columns must be understood as epoch theta1 theta2 dtheta1 dtheta2 action endpoint_y_position reward" << std::endl;
	  }

	};



	struct RBFParams {
	  static int nb_centers() { return 4;}

	  // The following defines the range for defining the centers of the basis functions
	  static std::pair<double, double> theta1_range() { return std::make_pair(-M_PI/4.0, 3.0 * M_PI/4.0);}
	  static std::pair<double, double> theta2_range() { return std::make_pair(-M_PI/4.0, 3.0 * M_PI_4);}
	  // The bounds for the speeds
	  static std::pair<double, double> dtheta1_range() { return std::make_pair(-9.0, 9.0);}
	  static std::pair<double, double> dtheta2_range() { return std::make_pair(-9.0, 9.0);}

	  // The followings define the ranges in which we switch to the LQR controller
	  static std::pair<double, double> theta1_range_for_lqr() { return std::make_pair(M_PI_2-M_PI/4.0, M_PI_2+M_PI/4.0);}
	  static std::pair<double, double> theta2_range_for_lqr() { return std::make_pair(-M_PI/4.0, M_PI/4.0);}
	  // The bounds for the speeds
	  static std::pair<double, double> dtheta1_range_for_lqr() { return std::make_pair(-M_PI_2, M_PI_2);}
	  static std::pair<double, double> dtheta2_range_for_lqr() { return std::make_pair(-M_PI_2, M_PI_2);}
	};

	double get_ith(int i, int Ni, std::pair<double, double> range) {
	  return range.first + double(i) * (range.second - range.first) / double(Ni-1.0);
	}

	template<typename STATE, typename PARAMS, typename AcrobatParams>
	struct Policy {

	  static int dimension() {
	    return PARAMS::nb_centers() * PARAMS::nb_centers() * PARAMS::nb_centers() * PARAMS::nb_centers()+1;
	  }
	  static double lbound(size_t index) { return -300;}
	  static double ubound(size_t index) { return  300;}

	  static double sample_rbf(STATE& s, double * params) {
	    double str = 0.0;
	    int param_index = 0;

	    double mu_theta1, mu_theta2, mu_dtheta1, mu_dtheta2;
	    double sigma_theta1 = 1.0 / double(PARAMS::nb_centers());
	    double sigma_theta2 = 1.0 / double(PARAMS::nb_centers());
	    double sigma_dtheta1 = (PARAMS::dtheta1_range().second -  PARAMS::dtheta1_range().first)/ double(2.0 * PARAMS::nb_centers());
	    double sigma_dtheta2 = (PARAMS::dtheta2_range().second -  PARAMS::dtheta2_range().first)/ double(2.0 * PARAMS::nb_centers());

	    str += params[param_index++];
	    double dist_theta1, dist_theta2, dist_dtheta1, dist_dtheta2;
	    unsigned int j_theta1, j_theta2, j_dtheta1, j_dtheta2;
	    for(j_theta1 = 0; j_theta1 < PARAMS::nb_centers(); ++j_theta1) {
	      mu_theta1 = get_ith(j_theta1, PARAMS::nb_centers(), PARAMS::theta1_range());
	      dist_theta1 = pow(sin(s._theta1 - mu_theta1), 2.0) / (2.0 * sigma_theta1 * sigma_theta1);
	      for(j_theta2 = 0; j_theta2 < PARAMS::nb_centers(); ++j_theta2) {
		mu_theta2 = get_ith(j_theta2, PARAMS::nb_centers(), PARAMS::theta2_range());
		dist_theta2 = pow(sin(s._theta2 - mu_theta2), 2.0) / (2.0 * sigma_theta2 * sigma_theta2);
		for(j_dtheta1 = 0; j_dtheta1 < PARAMS::nb_centers(); ++j_dtheta1) {
		  mu_dtheta1 = get_ith(j_dtheta1, PARAMS::nb_centers(), PARAMS::dtheta1_range());
		  dist_dtheta1 = pow(s._dtheta1 - mu_dtheta1, 2.0) / (2.0 * sigma_dtheta1 * sigma_dtheta1);
		  for(j_dtheta2 = 0; j_dtheta2 < PARAMS::nb_centers(); ++j_dtheta2) {
		    mu_dtheta2 = get_ith(j_dtheta2, PARAMS::nb_centers(), PARAMS::dtheta2_range());
		    dist_dtheta2 = pow(s._dtheta2 - mu_dtheta2, 2.0) / (2.0 * sigma_dtheta2 * sigma_dtheta2);
		    str += params[param_index++] * exp(- dist_theta1 
						       - dist_theta2
						       - dist_dtheta1
						       - dist_dtheta2);
		  }
		}
	      }
	    }
	    str =  AcrobatParams::max_torque() * tanh(str);
	    return str;
	  }

	  /**
	   * LQR controller with LQRD, TS=0.01, l2 = 1.0
	   **/
	  static double sample_lqr(STATE& s) {
	    double K[4] = {-240.7946 , -59.8844 , -88.3753 , -25.8847};
	    double str = -( K[0] * (s._theta1-M_PI_2) + K[1] * s._theta2 + K[2] * s._dtheta1 + K[3] * s._dtheta2);
	    bound(str, std::make_pair(-AcrobatParams::max_torque(), -AcrobatParams::max_torque()));
	    return str;
	  }
	  
	  static void sample_action(STATE& s, double* params, Action& a) {
	    double str;
	    if(is_in_bounds(s._theta1, PARAMS::theta1_range_for_lqr()) &&
	       is_in_bounds(s._theta2, PARAMS::theta2_range_for_lqr()) &&
	       is_in_bounds(s._dtheta1, PARAMS::dtheta1_range_for_lqr()) &&
	       is_in_bounds(s._dtheta2, PARAMS::dtheta2_range_for_lqr()))
	      a._f = sample_lqr(s);
	    else
	      a._f = sample_rbf(s, params);
	  }
	  
	};
	
      }
      
    }
  }
}


#endif
