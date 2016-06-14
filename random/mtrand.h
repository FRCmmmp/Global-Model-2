
#ifndef MTRAND_H
#define MTRAND_H

#include <vector>
#include <stdint.h>

class MTRand_int32 { // Mersenne Twister random number generator
public:
// default constructor: uses default seed only if this is the first instance
  MTRand_int32() { if (!init) seed(5489UL); init = true; }
// constructor with 32 bit int as seed
  MTRand_int32(uint32_t s) { seed(s); init = true; }
// constructor with array of size 32 bit ints as seed
  MTRand_int32(const uint32_t* array, int size) { seed(array, size); init = true; }
// the two seed functions
  void seed(uint32_t); // seed with 32 bit integer
  void seed(const uint32_t*, int size); // seed with array
// functions to set/get the state vector for checkpointing
  int get_state(std::vector<uint32_t>& iostate);
  void set_state(std::vector<uint32_t>& iostate, int32_t& iop);
// overload operator() to make this a generator (functor)
  uint32_t operator()() { return rand_int32(); }
// 2007-02-11: made the destructor virtual; thanks "double more" for pointing this out
  virtual ~MTRand_int32() {} // destructor
protected: // used by derived classes, otherwise not accessible; use the ()-operator
  uint32_t rand_int32(); // generate 32 bit random integer
private:
  static const int n = 624, m = 397; // compile time constants
// the variables below are static (no duplicates can exist)
  static uint32_t state[n]; // state vector array
  static int32_t p; // position in state array
  static bool init; // true if init function is called
// private functions used to generate the pseudo random numbers
  uint32_t twiddle(uint32_t, uint32_t); // used by gen_state()
  void gen_state(); // generate new state
// make copy constructor and assignment operator unavailable, they don't make sense
  MTRand_int32(const MTRand_int32&); // copy constructor not defined
  void operator=(const MTRand_int32&); // assignment operator not defined
};

// inline for speed, must therefore reside in header file
inline uint32_t MTRand_int32::twiddle(uint32_t u, uint32_t v) {
  return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)
    ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
}

inline uint32_t MTRand_int32::rand_int32() { // generate 32 bit random int
  if (p == n) gen_state(); // new state vector needed
// gen_state() is split off to be non-inline, because it is only called once
// in every 624 calls and otherwise irand() would become too big to get inlined
  uint32_t x = state[p++];
  x ^= (x >> 11);
  x ^= (x << 7) & 0x9D2C5680UL;
  x ^= (x << 15) & 0xEFC60000UL;
  return x ^ (x >> 18);
}

// generates double floating point numbers in the half-open interval [0, 1)
class MTRand : public MTRand_int32 {
public:
  MTRand() : MTRand_int32() {}
  MTRand(uint32_t seed) : MTRand_int32(seed) {}
  MTRand(const uint32_t* seed, int size) : MTRand_int32(seed, size) {}
  uint32_t i32() { return rand_int32(); }
  ~MTRand() {}
  double operator()() {
    return static_cast<double>(rand_int32()) * (1. / 4294967296.); } // divided by 2^32
private:
  MTRand(const MTRand&); // copy constructor not defined
  void operator=(const MTRand&); // assignment operator not defined
};

// generates double floating point numbers in the closed interval [0, 1]
class MTRand_closed : public MTRand_int32 {
public:
  MTRand_closed() : MTRand_int32() {}
  MTRand_closed(uint32_t seed) : MTRand_int32(seed) {}
  MTRand_closed(const uint32_t* seed, int size) : MTRand_int32(seed, size) {}
  ~MTRand_closed() {}
  double operator()() {
    return static_cast<double>(rand_int32()) * (1. / 4294967295.); } // divided by 2^32 - 1
private:
  MTRand_closed(const MTRand_closed&); // copy constructor not defined
  void operator=(const MTRand_closed&); // assignment operator not defined
};

// generates double floating point numbers in the open interval (0, 1)
class MTRand_open : public MTRand_int32 {
public:
  MTRand_open() : MTRand_int32() {}
  MTRand_open(uint32_t seed) : MTRand_int32(seed) {}
  MTRand_open(const uint32_t* seed, int size) : MTRand_int32(seed, size) {}
  ~MTRand_open() {}
  double operator()() {
    return (static_cast<double>(rand_int32()) + .5) * (1. / 4294967296.); } // divided by 2^32
private:
  MTRand_open(const MTRand_open&); // copy constructor not defined
  void operator=(const MTRand_open&); // assignment operator not defined
};

// generates 53 bit resolution doubles in the half-open interval [0, 1)
class MTRand53 : public MTRand_int32 {
public:
  MTRand53() : MTRand_int32() {}
  MTRand53(uint32_t seed) : MTRand_int32(seed) {}
  MTRand53(const uint32_t* seed, int size) : MTRand_int32(seed, size) {}
  ~MTRand53() {}
  double operator()() {
    return (static_cast<double>(rand_int32() >> 5) * 67108864. +
      static_cast<double>(rand_int32() >> 6)) * (1. / 9007199254740992.); }
private:
  MTRand53(const MTRand53&); // copy constructor not defined
  void operator=(const MTRand53&); // assignment operator not defined
};

#endif // MTRAND_H
