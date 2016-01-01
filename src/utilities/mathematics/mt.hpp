////////////////////////////////////////////////////////////////////////////////
//! \file		
//!     Mersenne Twister Algorithm
//!
//!		M. Matsumoto and T. Nishimura, "Mersenne Twister: A
//!		623-dimensionally equidistributed uniform pseudorandom number
//!		generator", ACM Trans. on Modeling and Computer Simulation Vol. 8,
//!		No. 1, January pp.3-30 (1998).
//!
//!		Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//!		All rights reserved.
//!
//!		Redistribution and use in source and binary forms, with or without
//!		modification, are permitted provided that the following conditions
//!		are met:
//!
//!		1. Redistributions of source code must retain the above copyright
//!		notice, this list of conditions and the following disclaimer
//!
//!		2. Redistributions in binary form must reproduce the above copyright
//!		notice, this list of conditions and the following disclaimer in the
//!		documentation and/or other materials provided with the distribution.
//!
//!		3. The names of its contributors may not be used to endorse or promote
//!		products derived from this software without specific prior written
//!		permission.	
//!
//!		THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//!		IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO
//!		, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//!		PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//!		CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//!		EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//!		PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//!		PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!		LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//!		NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//!		SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//!		http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//!	
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef _METRICS_MT_HPP_INCLUDED_
#define _METRICS_MT_HPP_INCLUDED_

////////////////////////////////////////////////////////////////////////////////
//!	\brief Mersenne Twister random number generator. See 
//!	http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html for details. 
//!
////////////////////////////////////////////////////////////////////////////////

class MersenneTwister
{
    public:
    
    MersenneTwister(void);
    ~MersenneTwister(void);

    double random(void) { return genrand_real1(); }
    void print(void);

    void init_genrand(unsigned long s);
    void init_by_array(unsigned long* init_key, int key_length);

    unsigned long genrand_int32(void);
    long genrand_int31(void);
    double genrand_real1(void);
    double genrand_real2(void);
    double genrand_real3(void);
    double genrand_res53(void);

    private:
    
    static const int N                    = 624;
    static const int M                    = 397;
    // constant vector a
    static const unsigned long MATRIX_A   = 0x9908b0dfUL;
    // most significant w-r bits
    static const unsigned long UPPER_MASK = 0x80000000UL;
    // least significant r bits
    static const unsigned long LOWER_MASK = 0x7fffffffUL;

    unsigned long* mt_;                  // the state vector
    int mti_;                            // mti == N+1 means mt not initialized

    unsigned long* init_key_;            // Storage for the seed vector
    int key_length_;                     // Seed vector length
    unsigned long s_;                    // Seed integer
    bool seeded_by_array_;               // Seeded by an array
    bool seeded_by_int_;                 // Seeded by an integer
};

#endif // METRICS_MT_H
