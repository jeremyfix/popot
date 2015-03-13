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

#ifndef POPOT_BENCHMARK
#define POPOT_BENCHMARK

#include <list>
#include <stdexcept>
#include <iterator>
#include<json/writer.h>

namespace popot
{
  namespace benchmark
  {
    /**
     * Benchmarking class; Given an algorithm and a number of trials for this algorithm
     * it computes the mean error, std, etc... and various statistics on the performances
     * We suppose ALGO provides :
     * - a constructor without arguments, 
     * - an init method
     * - a run method
     * - a step method, 
     * - a getBestFitness method,
     * - a fillBestPosition method
     * The Problem is also the one used as template for ALGO; It provides the stop method and the function 
     * evaluations counter
     *
     * Tips : To automatically estimate the number of trials
     * for getting a good approximate of the mean error of the algorithm, 
     * we may use the (Maurer & Pontil,2009) inequality with the sample mean and sample variance
     */

    template<typename ALGO, typename PROBLEM>
    class Benchmark
    {
    private:

      ALGO& _algo;
      PROBLEM& _problem;
      int _nb_trials;

      double _mean_error;
      double _std_error;
      int    _nb_failures;
      double * _best_position;
      double _best_fitness;
      double _success_rate;
      int _trial_counter;

      std::list<double> mean_error;
      std::list<double> std_error;
      std::list<double> mean_FE;


    public:
    Benchmark(ALGO& algo, PROBLEM& p, int nb_trials) 
      : _algo(algo), _problem(p), _nb_trials(nb_trials), _mean_error(0.0),_std_error(0.0), 
	_nb_failures(0), _best_position(0), _best_fitness(0), _success_rate(0.0), _trial_counter(0)
      {
	_best_position = new double[_problem.getDim()];
      }

      ~Benchmark()
      {
	delete[] _best_position;
      }

      /**
       * Returns the mean error over _nb_trials
       */
      double getMeanError(void) const
      {
	return _mean_error;
      }

      /**
       * Returns the std error over _nb_trials
       */
      double getStdError(void) const
      {
	return _std_error;
      }

      /**
       * Returns the success rate over _nb_trials
       */
      double getSuccessRate(void) const
      {
	return _success_rate;
      } 

      /**
       * Returns the best position over _nb_trials
       */
      void getBestPosition(double * pos) const
      {
	memcpy(pos, _best_position, _problem._dimension*sizeof(double));
      } 

      /**
       * Returns the best fitness over _nb_trials
       */
      double getBestFitness(void) const
      {
	return _best_fitness;
      }
 
      /**
       * Runs the benchmarks and computes the statistics
       */
      void run(int mode=0)
      {
	_mean_error = 0.0;
	_std_error = 0.0;
	_nb_failures = 0;
	_best_fitness = 0;
	_success_rate = 0.0;
	_trial_counter = 0;

	// Some temporary variables used to compute online the mean, std, .. 
	double sum_error = 0;
	double sum_square_errors=0;
	double init_fitness = 0;
        
	std::list<double> sum_error_vec;
	std::list<double> sum_error2_vec;
	mean_FE.clear();

	for(int i = 1 ; i <= _nb_trials ; ++i)
	  {
	    _problem.init();
	    _algo.init();
	    init_fitness = _algo.getBestFitness();

	    auto mean_FE_iter = mean_FE.begin(); // used only when i != 1
	    auto sum_error_iter = sum_error_vec.begin();
	    auto sum_error2_iter = sum_error2_vec.begin();

	    if(i == 1) {// The first run fills in the mean_FE 
	      mean_FE.push_back(_problem.getFE());
	      double best_fitness = _algo.getBestFitness();
	      sum_error_vec.push_back(best_fitness);
	      sum_error2_vec.push_back(best_fitness*best_fitness);
	    }
	    else {
	      if(*mean_FE_iter != _problem.getFE()) {
		std::ostringstream ostr;
		ostr.str("");
		ostr << "Cannot measure the mean error over several runs if each step of the algorithm does not produce the same function evaluation step, got " << _problem.getFE() << " while expecting " << *mean_FE_iter << std::endl;
		throw std::runtime_error(ostr.str());
	      }
	      ++mean_FE_iter;

	      double best_fitness = _algo.getBestFitness();
	      *(sum_error_iter++) += best_fitness;
	      *(sum_error2_iter++) += best_fitness * best_fitness;
	    }


	    while(_problem.getFE() < _problem.max_fe()) {
	      _algo.step();
	      if(i == 1) {// The first run fills in the mean_FE
		mean_FE.push_back(_problem.getFE());
		double best_fitness = _algo.getBestFitness();
		sum_error_vec.push_back(best_fitness);
		sum_error2_vec.push_back(best_fitness*best_fitness);
	      }
	      else {
		if(*mean_FE_iter != _problem.getFE()){
		  std::ostringstream ostr;
		  ostr.str("");
		  ostr << "Cannot measure the mean error over several runs if each step of the algorithm does not produce the same function evaluation step, got " << _problem.getFE() << " while expecting " << *mean_FE_iter << std::endl;
		  throw std::runtime_error(ostr.str());
		}
		++mean_FE_iter;

		double best_fitness = _algo.getBestFitness();
		*(sum_error_iter++) += best_fitness;
		*(sum_error2_iter++) += best_fitness * best_fitness;
	      }	
	    }


	    // Update the statistics
	    double best_fitness = _algo.getBestFitness();
	    sum_error += best_fitness;
	    sum_square_errors += best_fitness*best_fitness;

	    if(_problem.has_failed(best_fitness))
	      _nb_failures++;
	    if(i == 1 || best_fitness < _best_fitness)
	      {
		_best_fitness = best_fitness;
		_algo.fillBestPosition(_best_position);
	      }

	    ++_trial_counter;

	    // Eventually print the statistics
	    if(mode)
	      {
		_mean_error = sum_error/double(i);
		_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(i) * _mean_error*_mean_error)/double(i));
		_success_rate = 100.*double(i-_nb_failures)/double(i);

		print(std::cout);
	      }	   
	  }
		     
	// Before leaving, ensure that the statistics are up-to-date;
	_mean_error = sum_error/double(_nb_trials);
	if(_nb_trials > 1)
	  _std_error = sqrt(sum_square_errors/double(_nb_trials-1.0) - sum_error * sum_error/double(_nb_trials * (_nb_trials-1.0)));
	else {
	  _std_error = 0.0;
	  std::cerr << "I cannot compute a standard deviation with a single trial" << std::endl;
	}
	_success_rate = 100.*double(_nb_trials-_nb_failures)/double(_nb_trials);

	// Compute the mean error and variance accross the trials
	// The std is :  sigma(t) = sqrt{ 1/(N-1) \sum_i (e_i(t) - \mu(t))^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 2.0 * 1/(N-1) \mu(t) \sum_i e_i(t) + N/(N-1) \mu^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 2.0 * 1/(N(N-1)) (\sum_i e_i(t))^2 + N/(N-1) \mu^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 1/(N(N-1)) (\sum_i e_i(t))^2
	auto sum_error_iter = sum_error_vec.begin();
	auto sum_error2_iter = sum_error2_vec.begin();
	auto mean_error_iter = std::back_inserter(mean_error);
	auto std_error_iter = std::back_inserter(std_error);
	for(unsigned int i = 0 ; i < sum_error_vec.size(); ++i, ++sum_error_iter, ++sum_error2_iter, ++std_error_iter, ++mean_error_iter) {
	  *mean_error_iter = *sum_error_iter / double(_nb_trials);
	  *std_error_iter = sqrt(*sum_error2_iter / double(_nb_trials-1.0) - 1.0 / double(_nb_trials * (_nb_trials-1.0)) * (*sum_error_iter) * (*sum_error_iter)) ;
	}
      }

      /**
       * Print the statistics
       */
      void print(std::ostream & os) const
      {
	os << "Trial=" << _trial_counter << ";";
	os << "Error(mean)= " << getMeanError() << ";";
	os << "Error(std)= " << getStdError() << ";";
	os << "Best_min= " << getBestFitness() << ";";
	os << "Success_rate= " << getSuccessRate() << "%;";
	os << std::endl;
      }

      void dump_json(std::string filename, std::string problem_name, std::string algo_name) const {
	std::cout << "Dumping the benchmark results in " << filename << std::endl;
	
	Json::Value results;   
	Json::Value FE(Json::arrayValue);
	for(auto v: mean_FE)
	  FE.append(v);

	Json::Value mean(Json::arrayValue);
	for(auto v: mean_error)
	  mean.append(v);

	Json::Value std(Json::arrayValue);
	for(auto v: std_error)
	  std.append(v);
	
	results["problem"]["name"] = problem_name;
	results["problem"]["dimension"] = _problem.getDim();
	  
	results["algorithm"] = algo_name;
	results["nb_trials"] = _nb_trials;
	results["FE"] = FE;
	results["mean"] = mean;
	results["std"] = std;
	
	std::ofstream outfile(filename.c_str());
	outfile << results;
	outfile.close();
      }

      friend std::ostream & operator <<(std::ostream & os, const Benchmark &b)
	{
	  b.print(os);
	  return os;
	}
    };

    template<typename ALGO, typename PROBLEM>
      Benchmark<ALGO, PROBLEM>
      make_benchmark(ALGO& algo, PROBLEM& pb, int nb_trials)
      {
	return Benchmark<ALGO, PROBLEM>(algo, pb, nb_trials);
      };
  } // benchmark
} // popot

#endif
