#pragma once
#include <vector>
#include "Term.h"

/*
  A polynomial in R[x_1, x_2,..., x_N].
  Essentially an ordered list of terms.
  When a term is added to the list, the list is sorted and reduced as
  much as possible.  The list always has at least one element.

  Ordering is specified in constructor.
  Valid options are: "lex", "grlex", and "grevlex".
*/
template<unsigned int N>
class Polynomial
{
public:
  Polynomial(const char* monomialOrdering);

  Polynomial(const char* polynomial, const char* monomialOrdering);

  const char* getOrdering() const;

  Term<N> leadingTerm() const;

  void operator+=(const Term<N>& term);

  void operator-=(Term<N> term);

  void operator-=(const Polynomial<N>& polynomial);

  Polynomial<N> operator-(const Polynomial<N>& polynomial) const;

  Polynomial<N> operator*(const Term<N>& monomial) const;

  void printP() const;

  bool isZero() const;

private:
  std::vector<Term<N>> terms;
  const char* ordering;

  void sortTerms();

  void reduceTerms();

  void printM(const Term<N>& t) const;
};