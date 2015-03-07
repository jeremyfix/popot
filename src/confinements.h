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


#ifndef POPOT_CONFINEMENT_H
#define POPOT_CONFINEMENT_H

namespace popot
{
  namespace confinement
  {

    template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
    void confine(VECTOR& pos, VECTOR& vel, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
    {
      // In case the position is out of the bounds
      // we reset the velocities
      for(size_t i = 0 ; i < pos.size() ; ++i)
	{
	  if(pos[i] < lbound(i))
	    {
	      pos[i] = lbound(i);
	      vel[i] = 0;
	    }
	  else if(pos[i] > ubound(i))
	    {
	      pos[i] = ubound(i);
	      vel[i] = 0;
	    }
	}
    }

    template<typename VECTOR, typename LBOUND_FUNC, typename UBOUND_FUNC>
    void confine_spso2011(VECTOR& pos, VECTOR& vel, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
    {
      // In case the position is out of the bounds
      // we bound the position and revert the velocity halved
      for(size_t i = 0 ; i < pos.size() ; ++i)
	{
	  if(pos[i] < lbound(i))
	    {
	      pos[i] = lbound(i);
	      vel[i] *= -0.5;
	    }
	  else if(pos[i] > ubound(i))
	    {
	      pos[i] = ubound(i);
	      vel[i] *= -0.5;
	    }
	}
    }

  }
}

#endif
