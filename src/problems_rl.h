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

namespace popot {

  namespace problems {
    /**
     * Problems in reinforcement learning
     **/
    namespace rl {

      namespace mountain_car {

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
	    _theta = std::get<1>(PARAMS::theta_range()) + popot::math::uniform_random(0, 1) * (std::get<2>(PARAMS::theta_range()) - std::get<1>(PARAMS::theta_range()));
	    _dtheta = std::get<1>(PARAMS::dtheta_range()) + popot::math::uniform_random(0, 1) * (std::get<2>(PARAMS::dtheta_range()) - std::get<1>(PARAMS::dtheta_range()));
	  }
	};
       
	struct Action {
	  double _f;
	};

	template<typename TPOLICY, typename StateParams, typename PendulumParams>
	class Simulator
	{
	public:
	  static void transition() {
	  }

	  static double evaluate(double * params) {
	    return 0;
	  }
	};

      }

      namespace inverted_pendulum {

      }

      namespace acrobat {
      }

    }
  }
}


#endif
