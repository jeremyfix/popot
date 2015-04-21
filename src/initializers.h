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

#ifndef POPOT_INITIALIZER_H
#define POPOT_INITIALIZER_H

// For various initialization procedures and bounds handling, see
// Helwig, Wanka (2008) "Theoretical Analysis of Initial Particle Swarm Behavior"
// In Proceedings of the 10th International Conference on Parallel Problem Solving from Nature, G. Rudolph et al. (Eds.),LNCS 5199, pages 889–898

namespace popot
{
  namespace initializer
  {
    namespace position
    {
      template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
	void zero (VECTOR& p, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
	{
	  // fills in with 0
	  for(size_t i = 0 ; i < p.size() ; ++i)
	    p[i] = 0.0;
	};

      template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
      void uniform_random (VECTOR& p, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
      {
	// fills in with Uniform[min, max]
	for(size_t i = 0 ; i < p.size() ; ++i)
	  p[i] = popot::math::uniform_random(lbound(i),ubound(i));
      };
    } // namespace position

    namespace velocity
    {
      template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
	void zero(VECTOR& p, VECTOR& v, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
      {
	// fills in with 0
	for(size_t i = 0 ; i < v.size() ; ++i)
	  v[i] = 0.0;
      };

      template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
	void half_diff(VECTOR& p, VECTOR& v, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
      {
	// returns (U(min,max) - xi)/2.0
	for(size_t i = 0 ; i < v.size() ; ++i)
	  v[i] = 0.5*(popot::math::uniform_random(lbound(i), ubound(i))-p[i]);
      };

      template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
	void spso2011(VECTOR& p, VECTOR& v, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
      {
	// returns U(min,max) - xi
	for(size_t i = 0 ; i < v.size() ; ++i)
	  v[i] = popot::math::uniform_random(lbound(i), ubound(i))-p[i];
      };
    } // namespace velocity
  } // namespace initializer
} // namespace popot

#endif
