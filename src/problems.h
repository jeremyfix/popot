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

#ifndef POPOT_PROBLEMS_H
#define POPOT_PROBLEMS_H

namespace popot
{

  /**
   * Definition of some benchmarking problems
   */
  namespace problems
  {


    class Base {
    private:
      size_t _count;
      size_t _dimension;

    public:
      Base(size_t dimension): _count(0), _dimension(dimension) {}

      virtual ~Base(void) {}

      size_t getDim(void) const { return _dimension;}
      size_t getFE(void) const { return _count; }

      double evaluate(void * x) {
	++_count;
	return (*this)(x);
      }

      virtual void init(void) {
	_count = 0;
      }

      /**
       * The lower bound for the parameter index
       **/
      virtual double get_lbound(size_t index)           = 0;

      /**
       * The upper bound for the parameter index
       **/
      virtual double get_ubound(size_t index)           = 0;

      /**
       * Given a fitness and an epoch, should we stop to find for a solution
       **/
      virtual bool   stop(double fitness, size_t epoch) = 0;

      /**
       * Given a fitness, is it a failure or not ?
       **/
      virtual bool   has_failed(double fitness)         = 0;

      /**
       * What is the maximum number of Function evaluation allowed for this problem ?
       **/
      virtual size_t max_fe(void)                       = 0;

      /**
       * What is the lowest fitness we must obtain to consider a solution to be a success
       **/ 
      virtual double min_fitness(void)                  = 0;

      /**
       * What is the value of the function ?
       **/
      virtual double operator()(void * x)               = 0;
    };

    /**
     * N-dimensional Ackley function
     * @brief \f$ 20 (1 - \exp(-0.2 * \sqrt{\frac{1}{N} \sum_{i=1}^{N} x_i^2})) + \exp(0) - \exp(\frac{1}{N} \sum_{i=1}^{N} \cos(2\pi x_i) )\f$
     *  Bounds [-30; 30]
     */
    class Ackley : public Base
    {
    public:

      Ackley(int dimension) : Base(dimension){}

      double get_lbound(size_t index)
      {
	return -32.0;
      }

      double get_ubound(size_t index)
      {
	return 32.0;
      }

      bool stop(double fitness, size_t epoch)
      {
	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
	double * params = (double*) x;
	double fit = 0.0;
	double cos_x = 0.0;
	double sq_x = 0.0;
	size_t dimension = getDim();
	for(size_t i = 0 ; i < dimension ; ++i)
	  {
	    sq_x += pow(params[i],2.0);
	    cos_x += cos(2.0 * M_PI * params[i]);
	  }
	fit = 20.0 * (1.0 - exp(-0.2 * sqrt(1.0/double(dimension) * sq_x)))
	  + exp(1) - exp(1.0 / double(dimension) * cos_x);
	return fit;
      }
    };

    /**
     * N-dimensional Quadric function
     * @brief \f$ \sum_{i=1}^{N} (\sum_{j=1}^{i} x_i)^2 \f$
     *  Bounds [-100; 100]
     */
    class Quadric : public Base
    {
    public:

      Quadric(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) { return -100; }

      double get_ubound(size_t index) {	return 100;  }

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
    	double tmp_fit = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  {
    	    tmp_fit += params[i];
    	    fit += tmp_fit*tmp_fit;
    	  }
    	return fit;
      }
    };


    /**
     * N-dimensional Griewank function
     * @brief \f$ 1 + \frac{1}{4000} \sum_{i=1}^{N} x_i^2 - \prod_i cos(\frac{x_i}{\sqrt{i}})\f$
     *        Bounds [-600;600]
     */
    class Griewank : public Base
    {
    public:
      Griewank(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) { return -600; }

      double get_ubound(size_t index) {	return 600;  }

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
    	double cos_x = 1.0;
    	double sq_x = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  {
    	    sq_x += pow(params[i],2.0);
    	    cos_x *= cos(params[i]/sqrt(double(i+1)));
    	  }
    	fit = 1.0 + 1.0/4000.0 * sq_x - cos_x;
    	return fit;
      }
    };

    /**
     * N-dimensional Sphere function
     * @brief \f$ \sum_{i=1}^{N} x_i^2\f$
     * Bounds [-100;100]
     */
    class Sphere : public Base
    {
    public:

      Sphere(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) {	return -100; }

      double get_ubound(size_t index) {	return 100; }

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  fit += pow(params[i],2.0);
    	return fit;
      }
    };

    /**
     * N-dimensional Quartic function with uniform noise
     * @brief \f$ \sum_{i=1}^{N} (i x_i^4 + n_i)\f$
     * \f$ n_i \sim U(0,1)\f$ , Bounds [-1.28,1.28]
     */
    class QuarticNoise : Base
    { 
    public:

      QuarticNoise(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) {	return 1.28;}

      double get_ubound(size_t index) {	return -1.28;}

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  fit += double(i+1)*pow(params[i],4.0) + popot::math::uniform_random(0.0,1.0);
    	return fit;
      }
    };

    /**
     * N-dimensional Rastrigin function
     * @brief \f$ \sum_{i=1}^{N} (x_i^2 + 10 (1 - cos(2\pi x_i)))\f$
     * Bounds [-5.12,5.12]
     */
    class Rastrigin : public Base
    {
    public:

      Rastrigin(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) {	return -5.12; }

      double get_ubound(size_t index) { return 5.12;  }

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  fit += pow(params[i],2.0) + 10.0*(1.0 - cos(2.0*M_PI*params[i]));
    	return fit;
      }
    };

    /**
     * N-dimensional Rosenbrock banana function
     * @brief \f$ \sum_{i=1}^{N-1} (100 (y_{i+1} - y_i^2)^2 + (y_i - 1)^2)\f$
     *        with \f$ y_i = x_i + o_i - 1\f$
     * Bounds [-30,30]
     */
    class Rosenbrock : public Base
    {
    public:

      Rosenbrock(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) {	return -30; }

      double get_ubound(size_t index) {	return 30;}
	
      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 100;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
    	double y_i, y_i_1;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension-1 ; ++i)
    	  {
    	    y_i = params[i];
    	    y_i_1 = params[i+1];
    	    fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
    	  }
    	return fit;
      }
    };

    // class StaticRosenbrock
    // {
    // public:

    //   static double get_lbound(size_t index)
    //   {
    // 	return -30;
    //   }

    //   static double get_ubound(size_t index)
    //   {
    // 	return 30;
    //   }
	
    //   static bool stop(double fitness, size_t epoch)
    //   {
    // 	return (fitness <= 100) || (epoch >= 10000);
    //   }

    //   static double evaluate(double * x, size_t dimension)
    //   {
    // 	double fit = 0.0;
    // 	double y_i, y_i_1;
    // 	for(size_t i = 0 ; i < dimension-1 ; ++i)
    // 	  {
    // 	    y_i = x[i];
    // 	    y_i_1 = x[i+1];
    // 	    fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
    // 	  }
    // 	return fit;
    //   }
    // };


    /**
     * N-dimensional Schwefel1_2 function
     * @brief \f$ \sum_{i=1}^{N} (\sum_{j=1}^{i} x_i)^2\f$
     * Bounds [-500,500]
     */
    class Schwefel1_2 : public Base
    {
    public:

      Schwefel1_2(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) { return -500; }

      double get_ubound(size_t index) { return 500; }

      bool stop(double fitness, size_t epoch) {
    	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
    	double sum = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
    	  {
    	    sum += params[i];
    	    fit += pow(sum, 2.0);
    	  }
    	return fit;
      }
    };

    /**
     * N-dimensional Schwefel function
     * @brief \f$ 418.9829 N - \sum_{i=1}^{N} x_i \sin(\sqrt{|x_i|}) \f$
     * Bounds [-500,500]
     */
    class Schwefel : public Base
    {
    public:
      Schwefel(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) { return -500; }

      double get_ubound(size_t index) { return 500; }

      bool stop(double fitness, size_t epoch)
      {
	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }


      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
	size_t dimension = getDim();
    	double fit = 418.9829 * dimension;
    	for(size_t i = 0 ; i < dimension ; ++i)
	  fit += -params[i]*sin(sqrt(fabs(params[i])));
    	return fit;
      }
    };

    /**
     * N-dimensional Salomon function
     * @brief \f$ 1 + 0.1 \sqrt{\sum_{i=1}^{N} x_i^2} - cos(2\pi \sum_{i=1}^{N} x_i^2)\f$
     * Bounds [-600,600]
     */
    class Salomon : public Base
    {
    public:

      Salomon(int dimension) : Base(dimension) {}

      double get_lbound(size_t index) { return -600; }

      double get_ubound(size_t index) {	return 600; }

      bool stop(double fitness, size_t epoch)
      {
	return (fitness <= min_fitness()) || (getFE() >= max_fe());
      }

      bool has_failed(double fitness) {
	return fitness > min_fitness();
      }

      size_t max_fe(void) {
	return 10000*getDim();
      }

      double min_fitness(void) {
	return 1e-4;
      }

      double operator()(void * x)
      {
    	double * params = (double*) x;
    	double fit = 0.0;
    	double sum_sq = 0.0;
	size_t dimension = getDim();
    	for(size_t i = 0 ; i < dimension ; ++i)
	  sum_sq += pow(params[i],2.0);
    	fit = -cos(2.0 * M_PI * sum_sq) + 0.1 * sqrt(sum_sq) + 1.0;
    	return fit;
      }
    };

    // /**
    //  * 2-dimensional Dropwave function
    //  * @brief \f$ -\frac{1 + \cos(12 \sqrt{x_0^2 + x_1^2} + x_1^2)}{0.5 (x_0^2 + x_1^2) + 2} \f$
    //  * Bounds [-100,100]
    //  */
    // class Dropwave
    // {
    // public:
    //   size_t const dimension;
    //   size_t count;
      
    // Dropwave() : dimension(2), count(0) {}

    //   void init(void) {
    // 	count = 0;
    //   }
      
    //   double get_lbound(size_t index) 
    //   { 
    // 	return -5.12;
    //   }

    //   double get_ubound(size_t index) 
    //   { 
    // 	return 5.12;
    //   }

    //   bool stop(double fitness, size_t epoch) 
    //   {
    // 	return (fitness <= -1.0+1e-4) || (count >= 20000);
    //   }

    //   double evaluate(void * x){
    // 	double * params = (double*) x;
    // 	count++;
    // 	return -(1.0 + cos(12.0*sqrt(pow(params[0],2.0) + pow(params[1],2.0))))
    // 	  / (0.5 * (pow(params[0],2.0) + pow(params[1],2.0))+2.0);
    //   }
    // };


    // // The test functions of the CEC2005 benchmark
    // // "Problem Definitions and Evaluation Criteria for the CEC 2005 Special Session on Real-Parameter Optimization", 2005
    // // P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y. -P. Chen, A. Auger, S. Tiwari
    // namespace CEC2005
    // {
    // }

    // /**
    //  * Functions from BBOB : http://coco.lri.fr/downloads/download11.06/bbobc.tar.gz
    //  */
    // namespace BBOB
    // {
    // }

    // /**
    //  * Functions as defined within the Standard PSO 2011 (17/09/2012)
    //  */ 
    // namespace SPSO2011Bench
    // {

    //   /**
    //    * N-dimensional Rosenbrock banana function
    //    * @brief \f$ \sum_{i=1}^{N-1} (100 (y_{i+1} - y_i^2)^2 + (y_i - 1)^2)\f$
    //    *        with \f$ y_i = x_i + o_i - 1\f$
    //    * Bounds [-30,30]
    //    */
    // 	class Rosenbrock
    // 	{
    // 	public:
    // 	  size_t const dimension;
    // 	  size_t count;

    // 	Rosenbrock(int dimension) : dimension(dimension), count(0) {}

    // 	  void init(void) {
    // 	    count = 0;
    // 	  }

    // 	  double get_lbound(size_t index)
    // 	  {
    // 	    return -30;
    // 	  }

    // 	  double get_ubound(size_t index)
    // 	  {
    // 	    return 30;
    // 	  }
	
    // 	  bool stop(double fitness, size_t epoch)
    //       {
    // 	  return (fitness <= 1e-10) || (count >= 10000*dimension);
    //       }

    // 	  double evaluate(void * x)
    // 	  {
    // 	    double * params = (double*) x;
    // 	    count++;
    // 	    double fit = 0.0;
    // 	    double y_i, y_i_1;
    // 	    for(size_t i = 0 ; i < dimension-1 ; ++i)
    // 	      {
    // 		y_i = params[i]+1;
    // 		y_i_1 = params[i+1]+1;
    // 		fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
    // 	      }
    // 	    return fit;
    // 	  }
    // 	};

    //}

  } // namespace problems
} // namespace popot

#endif // POPOT_PROBLEMS_H
