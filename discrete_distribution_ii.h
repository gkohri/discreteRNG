/*
 *  
 *  Copyright (c) 2015, G.A. Kohring
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions 
 *  are met:
 *  
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *  
 *      * Redistributions in binary form must reproduce the above copyright 
 *        notice, this list of conditions and the following disclaimer in the 
 *        documentation and/or other materials provided with the distribution.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 *  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 *  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef DISCRETE_DISTRIBUTION_II_H
#define DISCRETE_DISTRIBUTION_II_H

/*
 * A class for producing random integers i in [0,n) with prbability p(i).
 *
 * This implementation uses the square histrogram method in combination with a 
 * small look-up table.
 *
 * (see: "Fast Generation of Discrete Random Variables", by G. Marsaglia, 
 *       W. W. Tsang and J. Wang, Journal of Statistical Software,
 *       July 2004, Vol. 11., Nr. 3. )
 */

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <map>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace rng {

template<typename _IntType = int >
class discrete_distribution_30bit {

    static_assert( std::is_integral<_IntType>::value,
                   "template argument not an integral type");


 private:

    std::vector<double> P;
    std::vector<double> V;
    std::vector<int>    K;
    int _J[256];

    int size;

 public:

    typedef _IntType result_type;

    /* 
     * The constructor takes as input iterators over probabilities. The 
     * probabilities should sum to one.
     *
     * Only the first thirty bits in the probabilities are used.
     */

    template <class InputIterator>
        discrete_distribution_30bit( InputIterator first, InputIterator last ): 
                                                          P( first, last ) {
   
        if ( P.size() == 0 ) {
            perror( "no probablities." );
            abort();
        }

        size = P.size();

        K.resize( size, 0 );
        V.resize( size, 0.0 );

        double max_val = static_cast<double>( ( (1 << 30) - 1 ) );
        double norm = 0.0;

        for ( int n = 0; n < size; ++n ) {
            P[n] = std::floor( ( P[n]*max_val ) );
            norm += P[n];
        }

        for ( int n = 0; n < size; ++n ) P[n] /= norm;

        init();
    }

    ~discrete_distribution_30bit(){}

    /* 
     * Generates the next random number in the sequence.
     *
     */
    template<typename _UniformRandomNumberGenerator>
    result_type operator()( _UniformRandomNumberGenerator& _urng ) {

        unsigned uran = _urng();

        int d = _J[ uran & 255 ];

        if ( d >= 0 ) {
            return d;
        } else {
            double U = static_cast<double>(uran)*2.328306437e-10;
            d = static_cast<int>( static_cast<double>( size )*U );
            if ( U < V[d] ) {
                return d;
            } else {
                return K[d];
            }
        }

    }


    void reset() {}

    /* 
     * Returns a vector containing the probabilities used in this
     * instance.
     *
     */
    std::vector<double> probabilities() const {
        return std::vector<double>( P.begin(), P.end() );
    }

    /* 
     * Returns the minimum random integer.
     *
     */
    result_type min() const { 
        return result_type( 0 );
    }

    /* 
     * Returns the maximum random integer.
     *
     */
    result_type max() const { 
        return result_type( size - 1 );
    }

    /* 
     * Check for equality by checking that the probabilities are the same.
     *
     */
    friend bool operator==( const discrete_distribution_30bit& v1, 
                            const discrete_distribution_30bit& v2) {
        return v1.P == v2.P;
    }


private:
    /* 
     * Initialize the look-up table and the histogram.
     *
     */
    void init() {

        std::vector<double> probs( P.begin(), P.end() );

        // ensure the probabilities are normalized
        double sum = 0.0;
        for ( int i = 0; i < size; ++i ) sum += probs[i];
        for ( int i = 0; i < size; ++i ) probs[i] /= sum;


        // fill in the look-up table
        int L = 0;
        double p_sum = 0.0;
        double  a = 1.0/static_cast<double>( size );
        for ( int i = 0; i < size; ++i ) {
            int k = static_cast<int>( 256.0*probs[i] );
            probs[i] = 256.0*probs[i] - static_cast<double>( k );
            p_sum += probs[i];
            for ( int j = 0; j < k; ++j ) _J[L+j] = i;
            L += k;
        }

        // fill empty table slots with -1
        for( int i = L; i < 256; ++i ) _J[i] = -1;

/*
        fprintf(stderr,"Table is %5.2f percent full\n",
                                    100.0*static_cast<double>(L)/256.0 );
*/


        // Initialize K, V
        for ( int j = 0; j < size; ++j ) {
            K[j] = j;
            V[j] = a*static_cast<double>(j + 1);
        }
        
        // normalize new probs

        double p_norm = 1.0/p_sum;
        for ( int i = 0; i < size; ++i ) probs[i] *= p_norm;

        std::multimap<double,int> mp;
        for ( int i = 0; i < size; ++i ) {
            mp.insert( std::pair<double,int>( probs[i], i ) );
        }

        // Apply Robin Hood rule

        for ( int i = 1; i <= size; ++i ) {
            std::multimap<double,int>::iterator imin = mp.begin();
            std::multimap<double,int>::iterator imax = std::prev( mp.end() );

            double min_p = (*imin).first;
            int min_i = (*imin).second;

            double max_p = (*imax).first;
            int max_i = (*imax).second;

            V[min_i] = min_p + static_cast<double>(min_i)*a;
            K[min_i] = max_i;

            mp.erase( std::prev( mp.end() ) );
            mp.erase( mp.begin() );

            max_p = (max_p + min_p) - a;

            mp.insert( std::pair<double,int>( a, min_i ) );
            mp.insert( std::pair<double,int>( max_p, max_i ) );

        }

    }

};

template<typename _IntType>
inline bool operator!=( const rng::discrete_distribution_30bit<_IntType>& v1, 
                        const rng::discrete_distribution_30bit<_IntType>& v2) {
        return !(v1.P == v2.P);
}

};  // namespace rng

#endif // DISCRETE_DISTRIBUTION_II_H
