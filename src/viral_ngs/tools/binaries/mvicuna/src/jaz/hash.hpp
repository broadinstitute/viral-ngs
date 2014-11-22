/***
 *  $Id$
 **
 *  File: hash.hpp
 *  Created: Apr 02, 2011
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2011-2012 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef HASH_HPP
#define HASH_HPP

#include <inttypes.h>
#include <string.h>


/** File: hash.hpp
 */
namespace jaz {

  /** Class: djb32
   *  Functor implementing 32bit hash DJB2 with XOR. It is hard to say if this is
   *  a good or a bad function! Some claim it is good, but some tests show
   *  it is not: http://www.team5150.com/~andrew/noncryptohashzoo/DJB.html.
   */
  class djb32 {
  public:
      /** Constructor: djb32
       *
       *  Parameter:
       *  seed - Default seed from http://en.wikipedia.org/wiki/2147483647.
       */
      explicit djb32(uint32_t seed = 2147483647) : seed_(seed) { }

      /** Function: operator()
       *  Performs hashing.
       */
      uint32_t operator()(const char* s, unsigned int l) const {
	  const unsigned char* S = (const unsigned char*)s;
	  uint32_t len = l;
	  uint32_t hash = 5381 + seed_ + len;

	  for (; len & ~1; len -= 2, S += 2) {
	      hash = ((((hash << 5) + hash) ^ S[0]) * 33) ^ S[1];
	  }

	  if (len & 1) hash = ((hash << 5) + hash) ^ S[0];

	  return hash ^ (hash >> 16);
      } // operator()

      /** Function: operator()
       *  Performs hashing.
       */
      uint32_t operator()(const std::string& s) const {
	  return this->operator()(s.c_str(), s.size());
      } // operator()

  private:
      uint32_t seed_;

  }; // class djb32


  /** Class: murmur264
   *  Functor implementing 64bit hash Murmur2.
   */
  class murmur264 {
  public:
      /** Constructor: murmur264
       */
      explicit murmur264(uint64_t seed = 2147483647) : seed_(seed) { }

      /** Function: operator()
       *  Performs hashing.
       */
      uint64_t operator()(const char* s, unsigned int l) const {
	  const unsigned char* S = (const unsigned char*)s;
	  uint64_t len = l;

	  const uint64_t m = 0xc6a4a7935bd1e995LL;
	  const uint8_t r = 47;

	  uint64_t h = len + seed_;

	  for (; len >= 8; len -= 8, S += 8) {
	      uint64_t k = (*(uint64_t*)S) * m;
	      k ^= k >> r;
	      k *= m;
	      h = (h * m) ^ k;
	  } // for

	  switch (len & 7) {
	    case 7: h ^= (uint64_t)S[6] << 48;
	    case 6: h ^= (uint64_t)S[5] << 40;
	    case 5: h ^= (uint64_t)S[4] << 32;
	    case 4: h ^= (uint64_t)S[3] << 24;
	    case 3: h ^= (uint64_t)S[2] << 16;
	    case 2: h ^= (uint64_t)S[1] << 8;
	    case 1: h ^= (uint64_t)S[0];
		h *= m;
	  } // switch

	  h ^= h >> r;
	  h *= m;
	  h ^= h >> r;

	  return h;
      } // operator()

      /** Function: operator()
       *  Performs hashing.
       */
      uint64_t operator()(const std::string& s) const {
	  return this->operator()(s.c_str(), s.size());
      } // operator()

  private:
      uint64_t seed_;

  }; // class murmur264


  /** Class: rabin64
   *  Functor implementing Rabin fingerprint. This code is based (to some extent)
   *  on Java implementations available in the Internet (with junk removed
   *  and several bug-fixes).
   */
  class rabin64 {
  public:
      // p_ represents: x^64 + x^4 + x^3 + x + 1 (which I hope is irreducible in Z2)

      /** Constructor: rabin64
       */
      rabin64() : p_(0x000000000000001BLL), p_deg_(64) { prv_init__(); }

      /** Function: operator()
       *  Performs hashing.
       */
      uint64_t operator()(const char* s, unsigned int l) const {
	  const unsigned char* S = (const unsigned char*)s;

	  uint64_t h = 0;
	  unsigned int pos = l % 8;

	  unsigned int i = 0;
	  if (pos != 0) for (; i < pos; ++i) h = (h << 8) ^ S[i];

	  while (i < l) {
	      h = tab32_[h & 0xFF] ^
		  tab40_[(h >> 8) & 0xFF] ^
		  tab48_[(h >> 16) & 0xFF] ^
		  tab56_[(h >> 24) & 0xFF] ^
		  tab64_[(h >> 32) & 0xFF] ^
		  tab72_[(h >> 40) & 0xFF] ^
		  tab80_[(h >> 48) & 0xFF] ^
		  tab88_[(h >> 56) & 0xFF] ^
		  ((uint64_t)S[i] << 56) ^
		  ((uint64_t)S[i + 1] << 48) ^
		  ((uint64_t)S[i + 2] << 40) ^
		  ((uint64_t)S[i + 3] << 32) ^
		  ((uint64_t)S[i + 4] << 24) ^
		  ((uint64_t)S[i + 5] << 16) ^
		  ((uint64_t)S[i + 6] << 8) ^
		  (uint64_t)S[i + 7];
	      i += 8;
	  } // while

	  return h;
      } // operator

      /** Function: operator()
       *  Performs hashing.
       */
      uint64_t operator()(const std::string& s) const {
	  return this->operator()(s.c_str(), s.size());
      } // operator()

  private:
      uint64_t p_;
      unsigned int p_deg_;

      void prv_init__() {
	  uint64_t xp = (uint64_t)1 << (p_deg_ - 1);
	  uint64_t* mods = new uint64_t[p_deg_];

	  mods[0] = p_;

	  for (unsigned int i = 1; i < p_deg_; i++) {
	      mods[i] = mods[i - 1] << 1;
	      if (mods[i - 1] & xp) mods[i] ^= p_;
	  }

	  memset(tab32_, 0, 256);
	  memset(tab40_, 0, 256);
	  memset(tab48_, 0, 256);
	  memset(tab56_, 0, 256);
	  memset(tab64_, 0, 256);
	  memset(tab72_, 0, 256);
	  memset(tab80_, 0, 256);
	  memset(tab88_, 0, 256);

	  for (unsigned int i = 0; i < 256; ++i) {
	      unsigned int c = i;
	      for (unsigned int j = 0; (j < 8) && (c > 0); ++j) {
		  if (c & 1) {
		      tab32_[i] ^= mods[j];
		      tab40_[i] ^= mods[j + 8];
		      tab48_[i] ^= mods[j + 16];
		      tab56_[i] ^= mods[j + 24];
		      tab64_[i] ^= mods[j + 32];
		      tab72_[i] ^= mods[j + 40];
		      tab80_[i] ^= mods[j + 48];
		      tab88_[i] ^= mods[j + 56];
		  } // if
		  c >>= 1;
	      } // for j
	  } // for i

	  delete[] mods;
      } // prv_init__

      uint64_t tab32_[256];
      uint64_t tab40_[256];
      uint64_t tab48_[256];
      uint64_t tab56_[256];
      uint64_t tab64_[256];
      uint64_t tab72_[256];
      uint64_t tab80_[256];
      uint64_t tab88_[256];

  }; // class rabin64

} // namespace jaz

#endif // HASH_HPP
