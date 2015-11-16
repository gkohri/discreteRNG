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

#ifndef _ZIPF_MANDELBROT_H
#define _ZIPF_MANDELBROT_H

/*
 * A class for producing random integers distributed according to the
 * Zipf-Mandelbrot probability mass function:
 *
 *           p(k;N,q,s) = 1/( H(N,q,s)*(k+q)^s )
 *
 * where,     
 *           H(N,q,s) = sum_(n=1)^(N) 1/(n+q)^s
 *
 * and, s>1, q >= 0, N > 1
 *
 * When q = 0 this becomes the mass function for Zipf's Law
 * When N -> infinity this becomes the Hurwitz Zeta mass function
 * When N -> infinity and q = 0, this becomes the Riemann Zeta mass function
 *
 */

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <tr1/cmath>
#include <vector>

namespace rng {


template< template<typename> class Discrete_Dist, typename Int_Type = uint32_t >
class ZipfMandelbrotDistribution {

    static_assert( std::is_integral<Int_Type>::value,
                   "template argument not an integral type");

 private:

    double s_;
    uint32_t q_;
    uint32_t N_;

    std::vector<double> probs_;

    Discrete_Dist<Int_Type>* dd_;

 public:

    typedef Int_Type result_type;

    /**
     * Creates a new ZipfMandelbrotDistribution given s, q, N
     * 
     *           p(k;N,q,s) = 1/( H(N,q,s)*(k+q)^s )
     *
     * where,     
     *           H(N,q,s) = sum_(n=1)^(N) 1/(n+q)^s
     *
     * and, s>1, q >= 0, N > 1
     *
     * Only s needs to be specified. The default for q is 0.
     * 
     * If N=0 (default) the maximum value of N realizeable with a 32 bit
     * unsigned integer given the values of q and s is calculated.
     */
    explicit ZipfMandelbrotDistribution( const double   s,
                                         const uint32_t q = 0, 
                                         const uint32_t N = 0 ) : s_(s), 
                                                                   q_(q), 
                                                                   N_(N),
                                                                   probs_(1) {

        if ( !(s_ > 1.0) ) {
            fprintf( stderr, "s (%g) must be greater than 1.0.\n", s_);
            abort();
        }

        if ( N_ < 1 ) {
            double dmax = 4294967295.0;
            double zeta_s = std::tr1::riemann_zeta( s_ );
            N_ = static_cast<uint32_t>( 
                                    std::pow( (dmax/zeta_s ), ( 1.0/s_ ) ) );
        }

        probs_.resize( N_, 0 );

        double p_sum = 0.0;
        for ( uint32_t k = 1; k < N_+1; ++k ) {
            double prob = 1.0/std::pow( static_cast<double>( k + q_ ), s_ );
            p_sum += prob;
            probs_[k-1] = prob;
        }

        double p_norm = 1.0/p_sum;
        for ( uint32_t i = 0; i < N_; ++i ) {
            probs_[i] *= p_norm;
        }

        dd_ = new Discrete_Dist<Int_Type>( probs_.begin(), probs_.end());

    }

    ~ZipfMandelbrotDistribution() {
        if ( dd_ != NULL ) {
            delete dd_;
            dd_ = NULL;
        }
    }

    template<typename UniformRandomNumberGenerator>
    result_type operator()( UniformRandomNumberGenerator& urng ) {
         return (*dd_)( urng ) + 1 + q_;
    }

    void reset() {}

    std::vector<double> probabilities() const {
        return dd_->probabilities();
    }

    result_type min() const { 
        return result_type( 1 );
    }

    result_type max() const { 
        return static_cast<result_type>( N_ );
    }

    friend bool operator==( const ZipfMandelbrotDistribution& v1, 
                            const ZipfMandelbrotDistribution& v2) {
        return v1.probs_ == v2.probs_;
    }

    friend bool operator!=( const ZipfMandelbrotDistribution& v1, 
                            const ZipfMandelbrotDistribution& v2) {
        return !( v1.probs_ == v2.probs_ );
    }


private:

};

};  // namespace rng

#endif // _ZIPF_MANDELBROT_H
