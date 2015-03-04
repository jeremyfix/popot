#ifndef POPOT_TOOLS_H
#define POPOT_TOOLS_H


// TODO :
// The mean and variance of the FIFO can be iteratively updated which is O(1) instead of O(N) , but this depends on whether we need to access frequently to the mean/variance

namespace popot
{
  namespace tools
  {

    template<int SIZE>
    class FIFO
    {
    private:
      class Element
      {
      private:
	Element * _next_element;
	Element * _previous_element;
	double _value;
      public:
	Element(double value)
	  {
	    _next_element = 0;
	    _previous_element = 0;
	    _value = value;
	  }

	Element * operator++(void)
	  {
	    return _next_element;
	  }

	double operator()(void)
	{
	  return _value;
	}

	double getValue(void)
	{
	  return _value;
	}

	void setValue(double value)
	{
	  _value = value;
	}

	bool hasNext(void)
	{
	  return _next_element != 0;
	}

	Element * getNext(void)
	{
	  return _next_element;
	}

	Element * getPrevious(void)
	{
	  return _previous_element;
	}

	void setNext(Element * elem)
	{
	  _next_element = elem;
	}

	void setPrevious(Element * elem)
	{
	  _previous_element = elem;
	}
	
      };

      Element * _first;
      Element * _last;
      int _size;
      
    public:
      FIFO(void)
      {
	_size = 0;
	_first = 0;
	_last = 0;
      }

      ~FIFO(void)
      {
	clear();
      }
      
      Element * begin(void)
      {
	return _first;
      }

      Element * end(void)
      {
	return 0;
      }

      void clear(void)
      {
	int nb_elem = _size;
	for(int i = 0 ; i < nb_elem ; ++i)
	  deleteFirst();
      }

      void insert(double value)
      {
	if(_size == 0)
	  {
	    _first = new Element(value);
	    _last = _first;
	    _size++;
	  }
	else if(SIZE == 1)
	  {
	    Element * first = _first;
	    _first = new Element(value);
	    _last = _first;
	    delete first;
	  }
	else if(_size < SIZE)
	  {
	    Element * elem = new Element(value);
	    _last->setNext(elem);
	    _last = elem;
	    _size++;
	  }
	else
	  {
	    Element * first = _first;
	    _first = _first->getNext();
	    Element * elem = new Element(value);
	    _last->setNext(elem);
	    _last = elem;
	    delete first;
	  }
      }

      void deleteFirst(void)
      {
	if(_size == 0)
	  throw popot::Exception::EmptyQueue();
	
	Element * first = _first;
	_first = _first->getNext();
	_size--;
	delete first;
      }

      int size(void)
      {
	return _size;
      }

      void print(void)
      {
	if(_size == 0)
	  {
	    std::cout << "Empty queue" << std::endl;
	  }
	else
	  {
	    Element * cur = _first;
	    std::cout << _size << " [";
	    for(int i = 0 ; i < _size ; ++i)
	      {
		std::cout << cur->getValue() << ";";
		cur = cur->getNext();
	      }

	    std::cout << "]" << std::endl;
	  }
      }

      double mean(void)
      {
	if(_size == 0)
	  throw popot::Exception::EmptyQueue();

	double mu = 0.0;
	Element * cur = _first;
	for(int i = 0 ; i < _size ; ++i)
	  {
	    mu += cur->getValue();
	    cur = cur->getNext();
	  }

	return mu / double(_size);
      }
      
      double variance(void)
      {
	if(_size == 0)
	  throw popot::Exception::EmptyQueue();

	double mu = mean();
	double var = 0.0;
	Element * cur;
	int i ;
	for(cur = _first, i = 0 ; i < _size ; ++i, cur=cur->getNext())
	  var += (cur->getValue() - mu) * (cur->getValue() - mu);

	var /= double(_size);
	return sqrt(var);
      }

      double min(void)
      {
	if(_size == 0)
	  throw popot::Exception::EmptyQueue();
	Element *iter;
	int i;
	double minval = _first->getValue();
	for(iter = begin() , i = 0; i < _size ; ++i, iter =iter->getNext())
	  minval = std::min(minval, (*iter)());
	return minval;
      }

      double max(void)
      {
	if(_size == 0)
	  throw popot::Exception::EmptyQueue();
	Element *iter;
	int i;
	double maxval = _first->getValue();
	for(iter = begin() , i = 0; i < _size ; ++i, iter =iter->getNext())
	  maxval = std::max(maxval, (*iter)());
	return maxval;
      }
      
    };

  }
}

#endif
