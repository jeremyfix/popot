#ifndef POPOT_EXCEPTIONS_H
#define POPOT_EXCEPTIONS_H

#include <sstream>

namespace popot
{
  class Exception
  {
  public :
    class Any {
    public :
      virtual ~Any(void) {}
      virtual std::string what(void) { return "";}
    };

    class EmptyQueue : public Any
    {
    public:
    EmptyQueue(void) : Any() {}
      virtual std::string what(void)
	{
	  return "Empty queue exception ";
	}
    };

    class MemoryNotAllocated: public Any
    {
    private:
      std::string _varname;
    public:
      MemoryNotAllocated(std::string varname) : Any() {
	_varname = varname;
      }
      
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Memory for " << _varname << " not allocated ";
	return s.str(); ;
      }
    };   

    class DimensionMismatch: public Any
    {
    private:
      int _d1, _d2;
    public:
      DimensionMismatch(int d1, int d2) : Any() {
	_d1 = d1;
	_d2 = d2;
      }
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Dimension mismatch exception; " << _d1 << " != " << _d2;
	return s.str();
      }
    };   

    class IndexOutOfRange: public Any
    {
    private:
      int _asked, _limit;
    public:
      IndexOutOfRange(int asked, int limit) : Any() {
	_asked = asked;
	_limit = limit;
      }
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Index out of range exception; Asked " << _asked << " ; limit : " << _limit;
	return s.str();
      }
    };

    class ParticleNotFound: public Any
    {

    public:
      ParticleNotFound() : Any() {}
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Particle not found exception ";
	return s.str();
      }
    };

    class FindBestFromEmptyNeighborhood: public Any
    {

    public:
      FindBestFromEmptyNeighborhood() : Any() {}
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Exception : cannot get the best particle of an empty neighborhood ";
	return s.str();
      }
    };

    class BestParticleNotInitialized: public Any
    {

    public:
      BestParticleNotInitialized() : Any() {}
      virtual std::string what(void)
      {
	std::ostringstream s;
	s << "Exception : best particle not initialized, call findBest before updateBest ! ";
	return s.str();
      }
    };



  }; // class Exception
} // namespace popot

#endif // POPOT_EXCEPTIONS_H
