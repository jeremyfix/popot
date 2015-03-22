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

      std::list<double> history_error;
      std::list<double> history_std;
      std::list<int> history_FE;
      int _record_every;


    public:
    Benchmark(ALGO& algo, PROBLEM& p, int nb_trials, int record_every=-1) 
      : _algo(algo), _problem(p), _nb_trials(nb_trials), _mean_error(0.0),_std_error(0.0), 
	_nb_failures(0), _best_position(0), _best_fitness(0), _success_rate(0.0), _trial_counter(0), _record_every(record_every) {
	_best_position = new double[_problem.getDim()];

	if(_record_every != -1) {
	  history_FE.push_back(0);
	  int FE = _record_every;	  
	  while(FE < p.max_fe()) {
	    history_FE.push_back(FE);
	    FE += _record_every;
	  }
	}
	  
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

	if(_record_every != -1) 
	  std::cout << "I will record the fitness every " << _record_every << " FE " << std::endl;

	// Some temporary variables used to compute online the mean, std, .. 
	double sum_error = 0;
	double sum_square_errors=0;
	double init_fitness = 0;
        
	std::list<double> sum_error_vec;
	std::list<double> sum_error2_vec;

	if(_record_every == -1)
	  history_FE.clear();
	// otherwise, the history_FE has been already filled with 
	// the FE at which the fitness must be recorded

	// max_fe_init allows to determine the first FE for which we have valid data
	// in case with use a uniform span of the [0, max_fe] range
	int max_fe_init = 0;

	for(_trial_counter = 1 ; _trial_counter <= _nb_trials ; ++_trial_counter) // When looping, the iterators are rewinded to the begining of the collections
	  {

	    // We initialize the iterators with which to fill the collections of errors
	    auto sum_error_iter = sum_error_vec.begin();
	    auto sum_error2_iter = sum_error2_vec.begin();
	    
	    // for history_FE_iter, we use a const iterator
	    // if necessary the first trial will fill the collection
	    // but otherwise we just read the collection
	    auto history_FE_iter =  history_FE.cbegin();

	    _problem.init();
	    _algo.init();
	    double current_fitness = _algo.getBestFitness();
	    int current_fe = _problem.getFE();


	    double last_fitness = current_fitness;
	    int last_fe = current_fe;

	    max_fe_init = std::max(max_fe_init, current_fe);
	    //std::cout << "Max fe init : " << max_fe_init << std::endl;


	    /////////
	    // We have to handle specifically how to record the error after initialization
	    // If _record_every == -1 , 
	    //    If this is the first trial, keep the FE at which to make the records
	    //    otherwise check we are at the same FE
	    // Otherwise 
	    //    We move the history_FE_iter pointer so that it points to a value
	    //    larger than or equal to the current_fe
	    if(_record_every == -1) {
	      if(_trial_counter == 1) {// The first run fills in the history_FE 
		history_FE.push_back(current_fe);
		// We expand the collection of errors
	        sum_error_vec.push_back(current_fitness);
		sum_error2_vec.push_back(current_fitness*current_fitness);
	      }
	      else {
		// we check if the current FE is coherent with the recorded one
		if(*history_FE_iter != _problem.getFE()) {
		  std::ostringstream ostr;
		  ostr.str("");
		  ostr << "Cannot measure the mean error over several runs if each step of the algorithm does not produce the same function evaluation step, got " << current_fe << " while expecting " << *history_FE_iter << "; you should call Benchmark by specifying a _record_every != -1 which leads to linearly interpolate the FE and errors on a regular span of [0, max_fe]." << std::endl;
		  throw std::runtime_error(ostr.str());
		}
		// if this is coherent, we don't have to record the FE

		++history_FE_iter;
		// and we sum the errors to the previously recorded ones
	        *(sum_error_iter++) += current_fitness;
		*(sum_error2_iter++) +=  current_fitness*current_fitness;
	      }
	    }
	    else {
	      // We move the FE_iter to that it points to a FE value
	      // above the current_FE
	      while((*history_FE_iter < current_fe) && history_FE_iter != history_FE.end()) {
		++history_FE_iter;

		if(_trial_counter == 1) { // the first run has to fill these first entries with a dummy fitness
		  sum_error_vec.push_back(0);
		  sum_error2_vec.push_back(0);
		}
		else {
		  ++sum_error_iter;
		  ++sum_error2_iter;
		}
	      }

	      // and there is no fitness to compute and store.

	    }


	    /////////////
	    // main loop
	    if(_record_every != -1) 
	    while(_problem.getFE() < _problem.max_fe()) {
	      _algo.step();
	      
	      current_fitness = _algo.getBestFitness();
	      current_fe = _problem.getFE();

	      if(_record_every == -1) {
		if(_trial_counter == 1) {// The first run fills in the history_FE 
		  history_FE.push_back(current_fe);
		  // We expand the collection of errors
		  sum_error_vec.push_back(current_fitness);
		  sum_error2_vec.push_back(current_fitness*current_fitness);
		}
		else {
		  // we check if the current FE is coherent with the recorded one
		  if(*history_FE_iter != _problem.getFE()) {
		    std::ostringstream ostr;
		    ostr.str("");
		    ostr << "Cannot measure the mean error over several runs if each step of the algorithm does not produce the same function evaluation step, got " << current_fe << " while expecting " << *history_FE_iter << "; you should call Benchmark by specifying a _record_every != -1 which leads to linearly interpolate the FE and errors on a regular span of [0, max_fe]." << std::endl;
		    throw std::runtime_error(ostr.str());
		  }
		  // if this is coherent, we don't have to record the FE

		  ++history_FE_iter;
		  // and we sum the errors to the previously recorded ones
		  *(sum_error_iter++) += current_fitness;
		  *(sum_error2_iter++) +=  current_fitness*current_fitness;
		}
	      }
	      else {

		// We need to interpolate the fitnesses for all the FE in history_FE
		// up to current_FE, using last_fitness and current_fitness

		// [last_fe, last_fitness] -> [current_fe, current_fitness]
		while(*history_FE_iter < current_fe && history_FE_iter != history_FE.end()) {
		  double interpolated_fitness = double(*history_FE_iter - last_fe)/(current_fe - last_fe) * (current_fitness - last_fitness)  + last_fitness;
		  if(_trial_counter == 1) {
		    sum_error_vec.push_back(interpolated_fitness);
		    sum_error2_vec.push_back(interpolated_fitness*interpolated_fitness);

		  }
		  else {
		    *sum_error_iter += interpolated_fitness;
		    ++sum_error_iter;
		    *sum_error2_iter +=  interpolated_fitness*interpolated_fitness;
		    ++sum_error2_iter;
	
		  }
		  ++history_FE_iter;
		}

	      }

	      last_fe = current_fe;
	      last_fitness = current_fitness;
	    }


	    // Update the statistics
	    current_fitness = _algo.getBestFitness();
	    current_fe = _problem.getFE();

	    if(_problem.has_failed(current_fitness))
	      _nb_failures++;
	    if(_trial_counter == 1 || current_fitness < _best_fitness)
	      {
		_best_fitness = current_fitness;
		_algo.fillBestPosition(_best_position);
	      }

	    // Eventually print the statistics
	    if(mode)
	      {
		_mean_error = sum_error/double(_trial_counter);
		_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(_trial_counter) * _mean_error*_mean_error)/double(_trial_counter));
		_success_rate = 100.*double(_trial_counter-_nb_failures)/double(_trial_counter);

		print(std::cout);
	      }	   
	  }

	// Adjust the trial_counter
	// when the loop finished it was as _nb_trials + 1
	--_trial_counter;
		     
	// Before leaving, ensure that the statistics are up-to-date;
	_mean_error = sum_error/double(_nb_trials);
	if(_nb_trials > 1)
	  _std_error = sqrt(sum_square_errors/double(_nb_trials-1.0) - sum_error * sum_error/double(_nb_trials * (_nb_trials-1.0)));
	else {
	  _std_error = 0.0;
	  std::cerr << "I cannot compute a standard deviation with a single trial" << std::endl;
	}
	_success_rate = 100.*double(_nb_trials-_nb_failures)/double(_nb_trials);


	if(_record_every != -1) {
	  // We must make some cleanups at the beginning of the collections of errors
	  // and drop all the elements for which FE < max_fe_init
	  // as not all the trials saved elements in this period
	  auto history_FE_iter = history_FE.begin();
	  while((history_FE_iter != history_FE.end()) && (*history_FE_iter < max_fe_init))
	    ++history_FE_iter;
	  int nb_elements_to_suppress = std::distance(history_FE.begin(), history_FE_iter) ;

	  auto history_FE_iter_end = history_FE.begin();
	  std::advance(history_FE_iter_end, nb_elements_to_suppress);
	  history_FE.erase(history_FE.begin(), history_FE_iter_end);

	  auto sum_error_iter_end = sum_error_vec.begin();
	  std::advance(sum_error_iter_end, nb_elements_to_suppress);
	  sum_error_vec.erase(sum_error_vec.begin(), sum_error_iter_end);

	  auto sum_error2_iter_end = sum_error2_vec.begin();
	  std::advance(sum_error2_iter_end, nb_elements_to_suppress);
	  sum_error_vec.erase(sum_error2_vec.begin(), sum_error2_iter_end);	  
	}



	// Compute the mean error and variance accross the trials
	// The std is :  sigma(t) = sqrt{ 1/(N-1) \sum_i (e_i(t) - \mu(t))^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 2.0 * 1/(N-1) \mu(t) \sum_i e_i(t) + N/(N-1) \mu^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 2.0 * 1/(N(N-1)) (\sum_i e_i(t))^2 + N/(N-1) \mu^2}
	//                        = sqrt{ (1/(N-1) \sum_i e_i(t)^2) - 1/(N(N-1)) (\sum_i e_i(t))^2

	auto sum_error_iter = sum_error_vec.begin();
	auto sum_error2_iter = sum_error2_vec.begin();
	auto history_error_iter = std::back_inserter(history_error);
	auto history_std_iter = std::back_inserter(history_std);
	for(unsigned int i = 0 ; i < sum_error_vec.size(); ++i, ++sum_error_iter, ++sum_error2_iter, ++history_std_iter, ++history_error_iter) {
	  *history_error_iter = *sum_error_iter / double(_nb_trials);
	  if(_nb_trials > 1)
	    *history_std_iter = sqrt(*sum_error2_iter / double(_nb_trials-1.0) - 1.0 / double(_nb_trials * (_nb_trials-1.0)) * (*sum_error_iter) * (*sum_error_iter)) ;
	  else
	    *history_std_iter = 0;
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
	for(auto v: history_FE)
	  FE.append(v);

	Json::Value mean(Json::arrayValue);
	for(auto v: history_error)
	  mean.append(v);

	Json::Value std(Json::arrayValue);
	for(auto v: history_std)
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
      make_benchmark(ALGO& algo, PROBLEM& pb, int nb_trials, int record_every=-1)
      {
	return Benchmark<ALGO, PROBLEM>(algo, pb, nb_trials, record_every);
      };
  } // benchmark
} // popot

#endif
