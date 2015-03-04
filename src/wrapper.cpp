/*   This file is part of popot
 *
 *   Copyright (C) 2011-2012
 *
 *   Author : Hadrien Glaude, Jeremy Fix
 *
 *   Contributor :
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *   Contact : Jeremy.Fix@gmail.com
 *
 */


#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;
#include "popot.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

struct LBOUND_FUNC {
    boost::python::object const lbound_fun;

    LBOUND_FUNC(boost::python::object const& lbound) : lbound_fun(lbound) {
    }

    double operator()(size_t index) const {
        lbound_fun.ptr();
        double r = boost::python::extract<double>(lbound_fun(index));
        return r;
    }
};

struct UBOUND_FUNC {
    boost::python::object const ubound_fun;

    UBOUND_FUNC(boost::python::object const& ubound) : ubound_fun(ubound) {
    }

    double operator()(size_t index) const {
        return boost::python::extract<double>(ubound_fun(index));
    }
};

struct STOP_CRITERIA {
    boost::python::object const stop_fun;

    STOP_CRITERIA(boost::python::object const& stop) : stop_fun(stop) {
    }

    double operator()(double fitness, size_t epoch) const {
        return boost::python::extract<double>(stop_fun(fitness, epoch));
    }
};

struct COST_FUNCTION {
    typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
    boost::python::object const cost_fun;
    size_t nbDim;

    COST_FUNCTION(boost::python::object const& cost, size_t dim) : cost_fun(cost), nbDim(dim) {
    }

    double operator()(TVector& pos) const {
        auto ptr = pos.getValuesPtr();
        boost::python::list parameters;
        for (size_t i = 0; i < nbDim; ++i) {
            parameters.append(ptr[i]);
        }
        return boost::python::extract<double>(cost_fun(parameters));
    }
};

class Wrapper_Stochastic_PSO_2006 {
private:
    typedef popot::PSO::particle::Particle<> ParticleSPSO;
    typedef popot::PSO::algorithm::Base<
    LBOUND_FUNC,
    UBOUND_FUNC,
    STOP_CRITERIA,
    COST_FUNCTION,
    void(*)(std::vector<ParticleSPSO >&,
            std::vector< typename ParticleSPSO::NeighborhoodType *> &,
            std::map< size_t, std::vector<size_t> > &),
    void(*)(ParticleSPSO&),
    void(*)(ParticleSPSO&),
    void(*)(ParticleSPSO&),
    void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
    void(*)(typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
    void(*)(typename ParticleSPSO::VECTOR_TYPE&, typename ParticleSPSO::VECTOR_TYPE&, const LBOUND_FUNC&, const UBOUND_FUNC&),
    ParticleSPSO> algo_type;
public:
    size_t const nbDim;
    LBOUND_FUNC lbound;
    UBOUND_FUNC ubound;
    STOP_CRITERIA stop;
    COST_FUNCTION cost;
    algo_type algo;

    Wrapper_Stochastic_PSO_2006(size_t dim, const boost::python::object& lbound_func, const boost::python::object& ubound_func, const boost::python::object& stop_func, const boost::python::object& cost_func) :
    nbDim(dim),
    lbound(LBOUND_FUNC(lbound_func)),
    ubound(UBOUND_FUNC(ubound_func)),
    stop(STOP_CRITERIA(stop_func)),
    cost(COST_FUNCTION(cost_func, nbDim)),
    algo(popot::algorithm::stochastic_spso2006(
    nbDim, 
    lbound, 
    ubound, 
    stop, 
    cost))
    {}
    
    double bestFitness() const {
        return algo.getBestFitness();
    }
    
    boost::python::list bestParticule() {
        auto v = algo.getBest().getPosition();
        boost::python::list parameters;
        for (size_t i = 0; i < nbDim; ++i) {
            parameters.append(v[i]);
        }
        return parameters;
    }
    
    size_t getEpoch() const {
        return algo.getEpoch();
    }
    
    double step() {
        return algo.step();
    }
    
    void run(int verbose) {
        algo.run(verbose);
    }
};

void seed(unsigned long seed) {
    RNG_GENERATOR::rng_srand(seed);
    RNG_GENERATOR::rng_warm_up();
}

BOOST_PYTHON_MODULE(libPyPopot) {
    using namespace boost::python;
    def("seed", seed, (boost::python::arg("seed")));
    class_<Wrapper_Stochastic_PSO_2006 > ("Stochastic_PSO_2006", init< size_t, boost::python::object, boost::python::object, boost::python::object, boost::python::object > ((boost::python::arg("nbDim"), boost::python::arg("lbound"), boost::python::arg("ubound"), boost::python::arg("stop"), boost::python::arg("cost"))))
        .def("bestFitness", &Wrapper_Stochastic_PSO_2006::bestFitness)
        .def("bestParticule", &Wrapper_Stochastic_PSO_2006::bestParticule)
        .def("getEpoch", &Wrapper_Stochastic_PSO_2006::getEpoch)
        .def("step", &Wrapper_Stochastic_PSO_2006::step)
        .def("run", &Wrapper_Stochastic_PSO_2006::run, (boost::python::arg("verbose")))
        ;
    
}
