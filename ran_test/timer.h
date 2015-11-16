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

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <cmath>
#include <ctime>

namespace util {

/**
 * A class for timing sections of code.  It measures the real time and cpu
 * time between calls.
 */
class Timer{
 public:

    /**
     * Creates a new timer and starts it running.
     */
    Timer() {
        realClockID = CLOCK_REALTIME;

        clock_getres( realClockID, &realResolution );

        realResolutionNS = static_cast<long>(realResolution.tv_sec)*1000000000l
                            + static_cast<long>(realResolution.tv_nsec);

        clock_gettime( realClockID, &realLast);

        realLastNS = static_cast<long>(realLast.tv_sec)*1000000000l +
                                    static_cast<long>(realLast.tv_nsec);

        if ( clock_getcpuclockid( 0, &cpuClockID ) == 0 ){
            clock_getres( cpuClockID, &cpuResolution );

            cpuResolutionNS = 
                    static_cast<long>(cpuResolution.tv_sec)*1000000000l +
                                    static_cast<long>(cpuResolution.tv_nsec);

            clock_gettime( cpuClockID, &cpuLast);

            cpuLastNS = static_cast<long>(cpuLast.tv_sec)*1000000000l +
                                    static_cast<long>(cpuLast.tv_nsec);
        } else {
            cpuClockID = 0;
        }
    };

    ~Timer(){}

    /**
     * Retrieves the resolution of this timer in nanoseconds
     */
    long getResolutionNS(){
        return realResolutionNS > cpuResolutionNS ?
                                        realResolutionNS : cpuResolutionNS ;
    }

    /**
     * Returns the time elapsed in seconds since the last time this
     * method was called or, if this is the first time it has been
     * called, since the timer was created.
     */
    void elapsed( double &realTime, double &cpuTime){
        timespec now;

        if ( cpuClockID != 0 ) {
            clock_gettime( cpuClockID, &now);
            long nowNS = static_cast<long>(now.tv_sec)*1000000000l +
                                    static_cast<long>(now.tv_nsec);
            long diff = nowNS - cpuLastNS;
            cpuLastNS = nowNS;
            cpuTime = static_cast<double>(diff)*1.0e-9;
        } else {
            cpuTime = 0.0;
        }

        clock_gettime( realClockID, &now);
        long nowNS = static_cast<long>(now.tv_sec)*1000000000l +
                                    static_cast<long>(now.tv_nsec);
        long diff = nowNS - realLastNS;
        realLastNS = nowNS;
        realTime = static_cast<double>(diff)*1.0e-9;
    }

 private:
    clockid_t realClockID;
    clockid_t cpuClockID;
    timespec realResolution;
    timespec cpuResolution;
    timespec realLast;
    timespec cpuLast;
    long realResolutionNS;
    long cpuResolutionNS;
    long realLastNS;
    long cpuLastNS;

    Timer(const Timer&) = delete;
    Timer& operator=(const Timer&) = delete;
};

}  // namespace util

#endif // END UTIL_TIMER_H
