#ifndef DISCRETE_DISTRIBUTION_H
#define DISCRETE_DISTRIBUTION_H

/*
 * A class for producing random integers i in [0,n) with prbability p(i).
 *
 * This implementation uses the alias method. It works with two arrays, 
 * one holding the probabilities and another holding the aliases.
 * Vose's algorithm is use for creating the arrays. It should
 * be O(log(n)) times faster than Walker's original algorithm.
 *
 * This class is a C++ adoption of a Java implementation by
 * Keith Schwarz (htiek@cs.stanford.edu).
 */

#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

namespace rng {

template<typename _IntType = int >
class discrete_distribution {

    static_assert( std::is_integral<_IntType>::value,
                   "template argument not an integral type");


 private:

    std::vector<double> weights;
    std::vector<std::pair<_IntType,double>> table;

    std::uniform_int_distribution<_IntType> uid;
    std::uniform_real_distribution<double> urd;


 public:

    typedef _IntType result_type;

    /* 
     * The constructor takes as input an iterator over positive weights.
     * The weights will be normalized to produce probabilities.
     *
     */

    template <class InputIterator>
    discrete_distribution( InputIterator first, InputIterator last ): 
                                                weights( first, last ),
                                                table( weights.size() ),
                                                uid( 0, weights.size()-1 ),
                                                urd( 0.0, 1.0 ) {
   
        if ( weights.size() == 0 ) {
            perror( "no probablities." );
            abort();
        }

        init();
    }

    template<typename _UniformRandomNumberGenerator>
    result_type operator()( _UniformRandomNumberGenerator &_urng ) {

        _IntType index = uid( _urng );

        return ( ( urd( _urng ) < table[index].second ) ? 
                                            index : table[index].first );
    }


    void reset() {}

    std::vector<double> probabilities() const {
        if ( weights.empty() ) {
            return std::vector<double>(1.0);
        } else {
            return std::vector<double>(weights.begin(),weights.end());
        }
    }

    result_type min() const { 
        return result_type(0);
    }

    result_type max() const { 
        return weights.empty() ? result_type(0) : 
                                 result_type( weights.size() - 1 );
    }

    friend bool operator==( const discrete_distribution & _v1, 
                            const discrete_distribution &_v2) {
        return _v1.weights == _v2.weights;
    }


private:
    void init() {

        double sum = std::accumulate( weights.begin(), weights.end(), 0.0 );

        if ( sum == 0.0 ) {
            perror( "discrete_distribution():weights sum to zero" );
            abort();
        }

        std::vector<double> probabilities( weights.size() );

        for (int i = 0; i < weights.size(); ++i) {
            weights[i] /= sum;
            probabilities[i] = weights[i];
        }

        double p_size = static_cast<double>( probabilities.size() );
        double average = 1.0/p_size;

        std::deque<int> small;
        std::deque<int> large;

        for (int i = 0; i < probabilities.size(); ++i) {
            if ( probabilities[i] >= average )
                large.push_back(i);
            else
                small.push_back(i);
        }

        while ( !small.empty() && !large.empty() ) {

            int less = small.back();
            small.pop_back();

            int more = large.back();
            large.pop_back();

            table[less].second = probabilities[less] * p_size;
            table[less].first = more;

            probabilities[more] += probabilities[less] - average;

            if ( probabilities[more] >= average ) {
                large.push_back( more );
            } else {
                small.push_back( more );
            }

        }

        while ( !small.empty() ) {
            table[ small.back() ].second = 1.0;
            small.pop_back();
        }

        while ( !large.empty() ) {
            table[ large.back() ].second = 1.0;
            large.pop_back();
        }

    }

};

template<typename _IntType>
inline bool operator!=( const rng::discrete_distribution<_IntType> & _v1, 
                        const rng::discrete_distribution<_IntType> &_v2) {
        return !(_v1.weights == _v2.weights);
}

};  // namespace rng

#endif // DISCRETE_DISTRIBUTION_H
