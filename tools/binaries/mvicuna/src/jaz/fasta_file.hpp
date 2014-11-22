/***
 *  $Id: fasta_file.hpp 120 2007-06-04 18:42:10Z zola $
 **
 *  File: fasta_file.hpp
 *  Created: Jul 23, 2006
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2006-2009 Jaroslaw Zola
 *  Please see attached LICENSE file.
 */

#ifndef BIO_IO_FASTA_FILE_HPP
#define BIO_IO_FASTA_FILE_HPP

#include <cstddef>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <utility>


/** Main bIO namespace.
 */
namespace bIO {

  /** This class provides interface to read input in the FASTA format.
   */
  class FASTA_input {
  public:
      /**
       */
      typedef std::pair<std::string, std::string> value_type;

      /** Constructs FASTA_input over the input stream @a in. Next it seeks
       *  the first sequence in the stream and reads its header (if possible).
       *  @param in is an input stream.
       */
      explicit FASTA_input(std::istream& in)
	  : owner_(false), input_(&in), val_(value_type("", "")), buf_("") {
	  pre_read_();
      } // FASTA_input

      /** Constructs FASTA_input for the file @a name. Next it seeks
       *  first sequence in the stream and reads its header (if possible).
       *  @param name is a name of the FASTA file to be read.
       */
      explicit FASTA_input(const std::string& name)
	  : owner_(true), val_(value_type("", "")), buf_("") {
	  input_ = new std::ifstream(name.c_str());
	  pre_read_();
      } // FASTA_input

      /**
       */
      virtual ~FASTA_input() { if (owner_ == true) delete input_; }


      /** Resets stream (tries to rewind to its beginning). Next it seeks
       *  first sequence in the stream and reads its header (if possible).
       *  @return true on success and false otherwise.
       */
      bool reset() {
	  input_->clear();
	  input_->seekg(0, std::ios_base::beg);
	  val_ = value_type("", "");
	  buf_ = "";
	  pre_read_();
	  return input_->good();
      } // reset

      /** @return state of the FASTA_input underlying stream.
       */
      operator bool() const { return input_->good(); }

      /** @return pair consisting of the last read sequence, that is
       *  containing sequence's header as the first argument, and
       *  sequence's chain as the second argument.
       */
      const value_type& operator*() const { return val_; }

      /** Tries to read FASTA sequence from the input.
       *  @return true on success and false otherwise.
       */
      bool operator++() {
	  if (input_->good() == false) return false;

	  if (buf_.empty() == false) {
	      unsigned int len = buf_.size() - 1;
	      if (buf_[len] == '\r') buf_.resize(len);
	  }

	  val_.first = (buf_.c_str() + 1);
	  val_.second = "";

	  do {
	      std::getline(*input_, buf_);
	      if ((buf_.empty() == false) && (buf_[0] != ';') && (buf_[0] != '>')) {
		  unsigned int len = buf_.size() - 1;
		  if (buf_[len] == '\r') buf_.resize(len);
		  val_.second += buf_;
	      }
	      if ((buf_.empty() == false) && (buf_[0] == '>')) break;
	  }
	  while (input_->good());

	  return true;
      } // operator++


  private:
      FASTA_input(const FASTA_input&);
      void operator=(const FASTA_input&);

      void pre_read_() {
	  while (input_->good()) {
	      std::getline(*input_, buf_);
	      if ((buf_.empty() == false) && (buf_[0] == '>')) break;
	  }

	  if (buf_.empty() == false) {
	      unsigned int len = buf_.size() - 1;
	      if (buf_[len] == '\r') buf_.resize(len);
	  }
      } // pre_read_

      bool owner_;
      std::istream* input_;

      value_type val_;
      std::string buf_;

  }; // class FASTA_input



  /** This class provides interface to write output in the FASTA format.
   */
  class FASTA_output {
  public:
      /**
       */
      typedef std::pair<std::string, std::string> value_type;

      /** Constructs FASTA_output over the output stream @a out.
       */
      explicit FASTA_output(std::ostream& out)
	  : owner_(false), output_(&out),
	    out_iter_(std::ostream_iterator<char>(*output_)),
	    val_(value_type("", "")) {
      } // FASTA_output

      /** Constructs FASTA_output over the file @a name.
       */
      explicit FASTA_output(const std::string& name)
	  : owner_(true), output_(new std::ofstream(name.c_str())),
	    out_iter_(std::ostream_iterator<char>(*output_)),
	    val_(value_type("", "")) {
      } // FASTA_output

      /**
       */
      virtual ~FASTA_output() { if (owner_ == true) delete output_; }


      /** @return state of the underlying stream.
       */
      operator bool() const { return output_->good(); }

      /** @return reference to the pair with sequence which will be written
       *  in the next call to operator++().
       */
      value_type& operator*() { return val_; }

      /** Writes currently stored sequence to the output, performing
       *  required formatting.
       *  @return true on success, false otherwise.
       */
      bool operator++() {
	  *(++out_iter_) = '>';
	  std::copy(val_.first.begin(), val_.first.end(), out_iter_);
	  *(++out_iter_) = '\n';

	  const int L = 79;

	  std::string::iterator iter = val_.second.begin();
	  std::size_t l = val_.second.size() / L;

	  for (std::size_t i = 0; i < l; ++i, iter += L) {
	      std::copy(iter, iter + L, out_iter_);
	      *(++out_iter_) = '\n';
	  }

	  std::copy(iter, val_.second.end(), out_iter_);
	  *(++out_iter_) = '\n';

	  return output_->good();
      } // operator


  private:
      FASTA_output(const FASTA_output&);
      void operator=(const FASTA_output&);

      bool owner_;
      std::ostream* output_;
      std::ostream_iterator<char> out_iter_;

      value_type val_;

  }; // class FASTA_output

} // namespace bIO

#endif // BIO_IO_FASTA_FILE_HPP
