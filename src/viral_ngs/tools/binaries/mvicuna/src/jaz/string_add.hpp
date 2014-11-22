/***
 *  $Id$
 **
 *  File: string_add.hpp
 *  Created: Jun 03, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_STRING_ADD_HPP
#define JAZ_STRING_ADD_HPP

#include <cstddef>
#include <iostream>
#include <functional>
#include <locale>
#include <sstream>
#include <string>


namespace jaz {

  /** Split a string into a list of strings.
   *  @param pat is a string separator.
   *  @param s is a string to split.
   *  @param out is a model of OutputIterator
   *  where the output should be stored.
   */
  template <typename charT, typename traits, typename Alloc, typename Iter>
  void split(charT pat,
	     const std::basic_string<charT, traits, Alloc>& s, Iter out) {
      unsigned int pos = 0;

      for (unsigned int i = 0; i < s.size(); ++i) {
	  if (s[i] == pat) {
	      if (i - pos > 0) {
		  *(out++) =
		      std::basic_string<charT, traits, Alloc>(s, pos, i - pos);
	      }
	      pos = i + 1;
	  }
      }

      *(out++) =
	  std::basic_string<charT, traits, Alloc>(s, pos, s.size() - pos);
  } // split


  /** Join the list [@a first, @a last) of strings into a single string.
   *  @param pat is a separator.
   *  @param init is a prefix for the joint string.
   */
  template <typename Iter, typename charT, typename traits, typename Alloc>
  std::basic_string<charT, traits, Alloc>
  join(charT pat, Iter first, Iter last, const std::basic_string<charT, traits, Alloc>& init) {
      std::basic_string<charT, traits, Alloc> s(init);
      if (s.empty() == true) s = *(first++);
      for (; first != last; ++first) s += pat + *first;
      return s;
  } // join

  template <typename Iter> inline std::string join(char pat, Iter first, Iter last) {
      return join(pat, first, last, std::string(""));
  } // join


  /** Remove spaces from a string.
   *  @param s a string to be processed.
   *  @param loc a locale to assess if a given char is space.
   *  @return a new string without spaces.
   */
  template <typename charT, typename traits, typename Alloc>
  std::basic_string<charT, traits, Alloc>
  remove_space(const std::basic_string<charT, traits, Alloc>& s,
	       const std::locale& loc = std::locale::classic()) {
      std::basic_string<charT, traits, Alloc> ns;
      for (std::size_t i = 0; i < s.size(); ++i) {
	  if (std::isspace(s[i], loc) == false) ns.push_back(s[i]);
      }
      return ns;
  } // remove_space


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_uc_compare :
	public std::binary_function<std::basic_string<charT, traits, Alloc>,
				    std::basic_string<charT, traits, Alloc>,
				    int> {
  public:
      explicit basic_uc_compare(const std::locale& loc = std::locale::classic())
	  : ct_(std::use_facet<std::ctype<char> >(loc)) { }

      int operator()(const std::basic_string<charT, traits, Alloc>& s1,
		     const std::basic_string<charT, traits, Alloc>& s2) const {

	  std::size_t l1 = s1.size();
	  std::size_t l2 = s2.size();

	  if (l1 < l2) return -1; else if (l2 < l1) return 1;

	  charT c1, c2;

	  for (std::size_t i = 0; i < l1; ++i) {
	      c1 = ct_.toupper(s1[i]);
	      c2 = ct_.toupper(s2[i]);
	      if (c1 < c2) return -1; else if (c2 < c1) return 1;
	  }

	  return 0;
      } // operator()

  private:
      const std::ctype<charT>& ct_;

  }; // class basic_uc_compare

  /** Functor to compare upper case version of strings.
   */
  typedef basic_uc_compare<char> uc_compare;


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_to_upper :
	public std::unary_function<std::basic_string<charT, traits, Alloc>,
				   std::basic_string<charT, traits, Alloc> > {
  public:
      explicit basic_to_upper(const std::locale& L = std::locale::classic())
	  : ct_(std::use_facet<std::ctype<char> >(L)), s_() { }

      const std::basic_string<charT, traits, Alloc>&
      operator()(const std::basic_string<charT, traits, Alloc>& s) const {
	  s_ = s;
	  for (std::size_t i = 0; i < s.size(); ++i) s_[i] = ct_.toupper(s[i]);
	  return s_;
      } // operator

  private:
      const std::ctype<charT>& ct_;
      mutable std::basic_string<charT, traits, Alloc> s_;

  }; // class basic_to_upper

  /** Functor to transform a string into its upper case form.
   */
  typedef basic_to_upper<char> to_upper;


  /** Simple approximate string matching.
   */
  template <typename charT, typename traits, typename Alloc>
  std::pair<long int, bool> approx_match(const std::basic_string<charT, traits, Alloc>& T,
					 const std::basic_string<charT, traits, Alloc>& P,
					 unsigned int mm) {
      unsigned int n = T.size();
      unsigned int m = P.size();

      long int pos = std::string::npos;
      bool mult = false;

      long int l = n - m + 1;
      unsigned int cmm = mm + 1;

      for (long int i = 0; i < l; ++i) {
	  unsigned int t = 0;
	  unsigned int j = 0;

	  for (; j < m; ++j) {
	      if (T[i + j] != P[j]) t++;
	      if (t == cmm) break;
	  }

	  if (j == m) {
	      if (t < cmm - 1) {
		  pos = i;
		  mult = false;
		  cmm = t + 1;
	      } else if (t == cmm - 1) {
		  if (pos == std::string::npos) pos = i;
		  else mult = true;
	      }
	  }
      } // for i

      return std::make_pair(pos, mult);
  } // approx_match


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_to_string {
  public:

      template <typename T>
      std::basic_string<charT, traits, Alloc> operator()(const T& t) {
	  os_.clear();
	  os_.str("");
	  os_ << t;
	  return os_.str();
      } // operator()

  private:
      std::basic_ostringstream<charT, traits, Alloc> os_;

  }; // class basic_to_string

  /** Functor to transform any serializable object to a string.
   */
  typedef basic_to_string<char> to_string;


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_string_to {
  public:

      template <typename T>
      bool operator()(const std::basic_string<charT, traits, Alloc>& s, T& t) {
	  is_.clear();
	  is_.str(s + '\n');
	  is_ >> t;
	  return is_.good();
      } // operator()

  private:
      std::basic_istringstream<charT, traits, Alloc> is_;

  }; // class basic_string_to

  /** Functor to transform a string into any unserializable object.
   */
  typedef basic_string_to<char> string_to;

} // namespace jaz

#endif // JAZ_STRING_ADD_HPP
