#ifndef POPOT_BENCHMARK
#define POPOT_BENCHMARK

namespace popot
{
  namespace benchmark
  {
    /**
     * Benchmarking class; Given an algorithm and a number of trials for this algorithm
     * it computes the mean error, std, etc... and various statistics on the performances
     * We suppose ALGO provides :
     * - a constructor without arguments, 
     * - a step method, 
     * - a getBestFitness method,
     * - a getBestPosition method
     * - a getEpoch() 
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
      double _mean_nb_f_evaluations;
      double _log_progress;
      int    _nb_failures;
      double * _best_position;
      double _best_fitness;
      double _success_rate;
      int _trial_counter;


    public:
    Benchmark(ALGO& algo, PROBLEM& p, int nb_trials) 
      : _algo(algo), _problem(p), _nb_trials(nb_trials), _mean_error(0.0),_std_error(0.0), _mean_nb_f_evaluations(0.0),
	_log_progress(0.0), _nb_failures(0), _best_position(0), _best_fitness(0), _success_rate(0.0), _trial_counter(0)
      {
	_best_position = new double[_problem._dimension];
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
       * Returns the mean number of function evaluations to meet the stop criteria over _nb_trials
       */
      double getMeanFEvaluations(void) const
      {
	return _mean_nb_f_evaluations;
      }

      /**
       * Returns the log progress over _nb_trials
       * NOT FUNCTIONAL YET
       */
      double getLogProgress(void) const
      {
	return _log_progress;
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
	_mean_nb_f_evaluations = 0.0;
	_log_progress = 0.0;
	_nb_failures = 0;
	_best_fitness = 0;
	_success_rate = 0.0;
	_trial_counter = 0;

	//ALGO* algo;
	// Some temporary variables used to compute online the mean, std, .. 
	double sum_error = 0;
	double sum_square_errors=0;
	double sum_fe = 0;
	double sum_log = 0;
	double init_fitness = 0;
	
	for(int i = 1 ; i <= _nb_trials ; ++i)
	  {
	    _problem.init();
	    _algo.init();

	    _trial_counter++;
	    
	    init_fitness = _algo.getBest().getFitness();
	    _algo.run();

	    // Update the statistics
	    sum_error += _algo.getBest().getFitness();
	    sum_square_errors += _algo.getBest().getFitness()*_algo.getBest().getFitness();
	    sum_fe += _problem._count;
	    sum_log += (log(_algo.getBest().getFitness()) - log(init_fitness))/double(_problem._count);
	    if(_problem.has_failed(_algo.getBest().getFitness()))
	      _nb_failures++;
	    if(i == 1 || _algo.getBest().getFitness() < _best_fitness)
	      {
		_best_fitness = _algo.getBest().getFitness();
		memcpy(_algo.getBest().getPosition().getValuesPtr(), _best_position, _problem._dimension*sizeof(double));
	      }
	      
	    // Eventually print the statistics
	    if(mode)
	      {
		_mean_error = sum_error/double(i);
		_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(i) * _mean_error*_mean_error)/double(i));
		_mean_nb_f_evaluations = sum_fe/double(i);
		_log_progress = sum_log/double(i);
		_success_rate = 100.*double(i-_nb_failures)/double(i);

		print(std::cout);
	      }
	  }
		     
	// Before leaving, ensure that the statistics are up-to-date;
	_mean_error = sum_error/double(_nb_trials);
	_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(_nb_trials) * _mean_error*_mean_error)/double(_nb_trials));
	_mean_nb_f_evaluations = sum_fe/double(_nb_trials);
	_log_progress = sum_log/double(_nb_trials);
	_success_rate = 100.*double(_nb_trials-_nb_failures)/double(_nb_trials);
      }

      /**
       * Print the statistics
       */
      void print(std::ostream & os) const
      {
	os << "Trial=" << _trial_counter << ";";
	os << "Error(mean)= " << getMeanError() << ";";
	os << "Error(std)= " << getStdError() << ";";
	os << "FE(mean)= " << getMeanFEvaluations() << ";";
	os << "Log_progress(mean)= " << getLogProgress() << ";";
	os << "Best_min= " << getBestFitness() << ";";
	os << "Success_rate= " << getSuccessRate() << "%;";
	os << std::endl;
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
