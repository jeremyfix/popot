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
