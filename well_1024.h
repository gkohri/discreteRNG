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

#ifndef RNG_WELL_1024_H
#define RNG_WELL_1024_H

#include <cstdint>
#include <algorithm>

namespace rng {


/**
 * An efficient, 32 bit random number generator 
 * proposed by Panneton, L'Ecuryer and Matsumoto. 
 * Well_1024 has a period of 2^1024 - 1 (10^308) and is one of a family of
 * similar generators.  For a detailed
 * description see: F. Panneton, P. L'Ecuryer and M. Matsumoto,
 * "Improved Long-Period Generators Based on Linear Recurrences Modulo 2",
 * ACM Transactions on Mathematical Software, 32, pp. 1-16 (2006).
 *
 * This implementation was evaluated using version 3.31.0 of the 
 * Dieharder suite of statistical tests for random number generators 
 * and successfully passed each test.
 *
 * The template parameter _uinit_type must be an unsigned integral type.
 * Note that since the underlying generator is a 32 bit generator, the
 * resulting random numbers will have at most 32 bits.
 *
 */
template<typename uint_type>
class well_engine {

 static_assert( std::is_unsigned<uint_type>::value, 
                  "uint_type must be unsigned integral type" );

 public:

    /** The type of the generated random value. **/
    typedef uint_type result_type;


    /** The default seed. **/
    static constexpr result_type default_seed = 69069u;


    /**
     * Constructs a Well random number generator using the specified seed.
     */
    explicit well_engine(result_type __seed = default_seed ) {
        seed( __seed );
    }


    /**
     * Reseeds the random number generator. This will reset the internal state
     * using the specified seed.
     */
    void seed( result_type __seed = default_seed ) {

        state_i = 0;

    // Initialize the internal state using a simple, poor quality rng.

        uint32_t s = static_cast<uint32_t>( __seed );
        for ( uint32_t j = 0; j < 32; j++) {
            s = 1103515245*s + 1013904243;
            state[j] = s;
        }

    // Move away from the initial state...
        discard(1000);

    };


    /**
     *  Return the minium number generated.
     */
    static constexpr result_type min() { return 0; }


    /**
     *  Return the maximum number generated.
     */
    static constexpr result_type max() { return 4294967295ul;}


    /**
     *  Discard the next z random numbers.
     */
    void discard( unsigned long long z) {
        for (; z != 0ull; --z ) (*this)();
    }


    /**
     *  Generate a random unsigned integer uniformly in the range [0,2^32-1]
     */
    inline
    result_type operator()() {
        uint32_t z0 = state[(state_i+31) & 31] ;
        uint32_t z1 = (state[state_i]) ^ m3_pos(8, state[(state_i + 3) & 31]);
        uint32_t z2 = m3_neg(19, state[(state_i + 24) & 31]) ^ 
                      m3_neg(14, state[(state_i + 10) & 31]);
        state[state_i] = z1 ^ z2;
        state[(state_i+31) & 31] = 
                m3_neg(11,z0) ^ m3_neg(7,z1) ^ m3_neg(13,z2) ;
        state_i = (state_i + 31) & 31;
        return state[state_i];
    }


    /**
     *  Check two well_engines for equality. The two are equal if they have
     *  the same internal state.
     */
    friend bool
    operator==( const well_engine& lhs, const well_engine& rhs ) { 
        return equal( lhs.state, rhs.state );
    }


 private:

    uint32_t state_i;
    uint32_t state[32];

    inline uint32_t m3_neg(const uint32_t& t,const uint32_t& v) {
        return (v^(v<<(t)));
    }

    inline uint32_t m3_pos(const uint32_t& t,const uint32_t& v) {
        return (v^(v>>t));
    }

};


/**
 *  Check two well_engines for inequality.
 */
template<typename uint_type>
inline bool
operator!=( const rng::well_engine<uint_type>& lhs,
            const rng::well_engine<uint_type>& rhs ) { 
    return !( lhs == rhs );
};


/**
 *  The default RNG returns 32 bit unsigned integers
 */
typedef well_engine<uint32_t> Well_1024;


}  // namespace rng

#endif   // END RNG_WELL_1024_H
