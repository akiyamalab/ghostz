/*
 * alphabet_type.h
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#ifndef ALPHABET_TYPE_H_
#define ALPHABET_TYPE_H_

#ifdef __GNUC__
#include <tr1/memory>
#else
#include <string>
// for nvcc 
#define __aligned__ ignored
#include <boost/tr1/memory.hpp>
#undef __aligned__
#endif


class AlphabetType {
public:
  virtual ~AlphabetType();
  virtual std::string GetRegularLetters() const = 0;
  virtual std::string GetAmbiguousLetters() const = 0;
  virtual char GetUnknownLetter() const = 0;

};

typedef std::tr1::shared_ptr<AlphabetType> AlphabetTypePtr;

#endif /* ALPHABET_TYPE_H_ */
