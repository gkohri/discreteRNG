# discreteRNG

###A C++ library for generating discrete random numbers according to the Zipf-Mandelbrot distrbution.

## Zipf-Mandelbrot

The Zipf-Mandelbrot distribution finds uses in many applications. The general form of the probability mass function is:


> f(k;N,q,s) = 1/( H(N,q,s)*(k+q)^s )

where

> H(N,q,s) = sum_(n=1)^(N) 1/(n+q)^s

and, 

>  s,N > 1 and q >= 0


When q = 0 this becomes the mass function for Zipf's Law

When N -> infinity this becomes the Hurwitz Zeta mass function

When N -> infinity and q = 0, this becomes the Riemann Zeta mass function


## Algorithm

The approach taken here is to generate a table of probabilities according to the paramaters N, q and s, and then use a general purpose algorithm for generating discrete random numbers.

All classes  are designed to work within the C++11 standard framework for Pseudo-random number generation.

## Classes

`zipf_mandelbrot_distribution`: A template class for creating the random numbers. It takes as input, the values for N, q and s, and, as template parameter, a class for creating discrete random numbers. Although this class can work with the standard `discrete_distribution` template class from C++11, we provide two other templates that are more efficient.

`discrete_distribution`:  A class using Walker's histogram algorithm to generate discrete random numbers. The tables are created using Vose's algorithm.

`discrete_distribution_30bit`: This class uses the MTW (Marsaglia, Tsang and Wang) algorithm for creating the discrete random numbers. As the name suggests, the probabilities are limited to 30 bit accuracy.

`discrete_distribution_30bit` is about four times faster than `discrete_distribution` when generating random deviates, but needs much longer for the initiation. Both classes are faster than the standard implementation.

`well_1024`: A uniform, 32 bit, psuedo random number generator with a period of 2^1024. The Well algorith will thermalize faster than Matsumoto and Nishimura's Mercene twister. The implementation here is faster than the standard `mersenne_twister`.

## Tests and example usage

The directory **ran_test** contains a test program which demonstrates how to use the classes together with the C++11 framework for psuedo random number generation.

