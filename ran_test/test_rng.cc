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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iterator>
#include <random>
#include <vector>

#include "well_1024.h"
#include "discrete_distribution.h"
#include "discrete_distribution_ii.h"
#include "timer.h"
#include "zipf-mandelbrot.h"

#include "kfunc.h"

using std::vector;
using std::uniform_real_distribution;
using std::mt19937;

using rng::well_1024;
using rng::discrete_distribution_30bit;
using rng::zipf_mandelbrot_distribution;
using util::Timer;


template<class _DiscDist, class _RNG>
void test_disc_dist( _DiscDist &dd, _RNG &rng, const int num_gen ) {


    Timer timer;
    double real_time, cpu_time;

    vector<double> probs = dd.probabilities();

    vector<long> obsv( probs.size(), 0 );

    int offset = dd.min();

    // Generate the random numbers

    timer.elapsed( real_time, cpu_time );
    for ( long it = 0; it < num_gen; ++it ) {
        int index = dd( rng ) - offset;
        obsv[index] += 1;
    }
    timer.elapsed( real_time, cpu_time );


    // rate

    double gen_rate = static_cast<double>(num_gen)/real_time;
    fprintf( stdout,"Discrete integer generation: %14.7e rands/s \n",
                                            gen_rate );


    // statistics

    vector<double> expected( probs.size(), 0 );
    vector<bool> low( probs.size(), false );

    // ... expected number
    double d_num_gen = static_cast<double>( num_gen );
    for ( size_t n = 0; n < expected.size(); ++n ) {
        expected[n] = probs[n]*d_num_gen;
        if ( expected[n] < 10.0 ) low[n] = true;
    }

    double sum_dev = 0.0;
    double g_stat = 0.0;
    double low_bin_obs  = 0.0;
    double low_bin_exp  = 0.0;
    double low_bin_prob = 0.0;
    int df = 0;

    // ... good counts

    for ( int n = 0; n < probs.size(); ++n ) {
        if ( low[n] ) {
            low_bin_prob += probs[n];
            low_bin_obs += static_cast<double>( obsv[n] );
        } else {
            double expect = expected[n];
            double observ = static_cast<double>( obsv[n] );
            double dev = observ - expect;
            sum_dev += (dev*dev)/expect;
            ++df;
            if ( observ > 0.0 ) {
                g_stat += observ*std::log( (observ/expect) );
            }
        }
    }

    // ... low counts

    if ( low_bin_prob > 0.0 ) {
        low_bin_exp = low_bin_prob*d_num_gen;
        double dev = low_bin_obs - low_bin_exp;
        sum_dev += (dev*dev)/low_bin_exp;
        ++df;
        if ( low_bin_obs > 0.0 ) {
            g_stat += low_bin_obs*std::log( (low_bin_obs/low_bin_exp) );
        }
    }

    // ... correct likelihoods

    g_stat *= 2.0;

    // output results

    fprintf( stdout, "degress of freedom: %d\n", df ); 
    fprintf( stdout, "mean squared deviation: %14.7e\n", sum_dev ); 
    fprintf( stdout, "log likelihood: %14.7e\n", g_stat ); 
    fprintf( stdout, "p value msd: %14.7e\n", 
                kf_gammaq( ( std::floor(static_cast<double>(df-1)/2.0 )), 
                           ( sum_dev/2.0 ) ) ); 
    fprintf( stdout, "p value log likelihood: %14.7e \n", 
                kf_gammaq( ( static_cast<double>(df-1)/2.0 ),
                                        ( std::abs(g_stat)/2.0 ) ) ); 

/*
 * Uncomment the following code to see the individual probabilities
 *
*/

/*
    for ( int n = 0; n < obsv.size(); ++n ) {
        if ( !low[n] ) {
            fprintf( stdout, "%10d\t%10ld\t%15.1f\t%17.10e\n", n, obsv[n], 
                            expected[n], probs[n]);
        }
    }
    fprintf( stdout, "%10d\t%lf\t%lf\t%17.10e\n", df, low_bin_obs,
                            low_bin_exp, low_bin_prob );
*/

}

int main( int argc, char* argv[] ) {

    long num_gen = 100000000;
    int N = 1000;

    // choose a uniform random number generator

    well_1024 urng{ 15791113 };
    //mt19937 urng{ 15791113 };

    vector<double> weights(N, 0.0);

    fprintf( stdout,"Test Random table............. \n");

    // Set up a small table of random integers

    std::uniform_int_distribution<> uint( N, ( ( 1 << 30 ) - 1 ) );

    double w_sum = 0.0;
    for ( int n = 0; n < N; ++n ) {
        weights[n] = uint( urng );
        w_sum += weights[n];
    }
    
    double w_norm = 1.0/w_sum;
    for ( int n = 0; n < N; ++n ) {
        weights[n] *= w_norm;
    }

    Timer timer;

    double real_time, cpu_time;

    timer.elapsed( real_time, cpu_time );

    discrete_distribution_30bit<int> disc_dist(weights.begin(),weights.end());

    timer.elapsed( real_time, cpu_time );

    fprintf( stdout,"Framework generation: %14.7e s \n", real_time);

    test_disc_dist( disc_dist, urng, num_gen );

    // 2nd test: Zipf Mandelbrot

    fprintf( stdout,"\nTest Zipf Mandelbrot Distribution ............. \n");

    timer.elapsed( real_time, cpu_time );

    double alpha = 1.75;
    zipf_mandelbrot_distribution<discrete_distribution_30bit,int> zm( alpha );

    timer.elapsed( real_time, cpu_time );

    fprintf( stdout,"Framework generation: %14.7e s \n", real_time);


    test_disc_dist( zm, urng, num_gen );

}
