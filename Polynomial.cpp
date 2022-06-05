#include <iostream>
#include <cmath>
#include "Polynomial.h"

/*
  Creates 0 polynomial with specific ordering.
  Valid options are "lex", "grlex", and "grevlex".
*/
template<unsigned int N>
Polynomial<N>::Polynomial(const char* monomialOrdering) : ordering(monomialOrdering)
{
  terms.emplace_back(Term<N>());
}

/*
  Creates Polynomial object from a string.  Example:
  Polynomial("lex", "-x^3y^2z + 3z^3 - xy + z");
*/
template<unsigned int N>
Polynomial<N>::Polynomial(const char* polynomial, const char* monomialOrdering) : ordering(monomialOrdering)
{
  for (int i = 0; polynomial[i] != '\0';)
  {
    bool detectedTerm = false;
    bool minus = false;
    int coefficient = 1;
    std::array<unsigned int, N> alpha{};

    // Increments forward until start of new term is detected
    while (polynomial[i] != '-' && (polynomial[i] < 'z' - N + 1 || polynomial[i] > 'z') && (polynomial[i] < '1' || polynomial[i] > '9') && polynomial[i] != '\0') ++i;

    // Stores minus sign on coefficient
    if (polynomial[i] == '-')
    {
      minus = true;
      ++i;
    }

    // Stores coefficient and powers on x,y,z
    while (polynomial[i] != '+' && polynomial[i] != '-' && polynomial[i] != '\0')
    {
      if (polynomial[i] >= '1' && polynomial[i] <= '9')
      {
        coefficient = 0;
        detectedTerm = true;
        std::vector<int> digits = { polynomial[i] - '0' };
        ++i;
        while (polynomial[i] >= '0' && polynomial[i] <= '9')
        {
          digits.emplace_back(polynomial[i] - '0');
          ++i;
        }
        unsigned int K = digits.size();
        for (unsigned int i = 0; i < K; ++i)
          coefficient += digits[K - i - 1] * pow(10, i);
      }
      else if (polynomial[i] >= 'z' - N + 1 && polynomial[i] <= 'z')
      {
        detectedTerm = true;
        char variable = polynomial[i];
        ++i;
        if (polynomial[i] != '^') alpha[variable - ('z' - N + 1)] = 1;
        else
        {
          ++i;
          if (polynomial[i] >= '1' && polynomial[i] <= '9')
          {
            int exponent = 0;
            std::vector<int> digits = { polynomial[i] - '0' };
            ++i;
            while (polynomial[i] >= '0' && polynomial[i] <= '9')
            {
              digits.emplace_back(polynomial[i] - '0');
              ++i;
            }
            unsigned int K = digits.size();
            for (unsigned int i = 0; i < K; ++i)
              exponent += digits[K - i - 1] * pow(10, i);
            alpha[variable - ('z' - N + 1)] = exponent;
          }
        }
      }
      else ++i;
    }

    // Adds term to polynomial
    if (detectedTerm)
    {
      if (minus) coefficient *= -1;
      *this += Term<N>(coefficient, alpha);
    }
  }
}

template<unsigned int N>
const char* Polynomial<N>::getOrdering() const
{
  return ordering;
}

template<unsigned int N>
Term<N> Polynomial<N>::leadingTerm() const
{
  return terms[0];
}

template<unsigned int N>
void Polynomial<N>::operator+=(const Term<N>& term)
{
  terms.emplace_back(term);
  sortTerms();
  reduceTerms();
}

template<unsigned int N>
void Polynomial<N>::operator-=(Term<N> term)
{
  term.c *= -1;
  terms.emplace_back(term);
  sortTerms();
  reduceTerms();
}

template<unsigned int N>
void Polynomial<N>::operator-=(const Polynomial<N>& polynomial)
{
  for (unsigned int i = 0; i < polynomial.terms.size(); ++i)
    *this -= polynomial.terms[i];
}

template<unsigned int N>
Polynomial<N> Polynomial<N>::operator-(const Polynomial<N>& polynomial) const
{
  Polynomial p = *this;
  for (unsigned int i = 0; i < polynomial.terms.size(); ++i)
    p -= polynomial.terms[i];
  return p;
}

template<unsigned int N>
Polynomial<N> Polynomial<N>::operator*(const Term<N>& monomial) const
{
  Polynomial p = Polynomial(ordering);
  for (unsigned int i = 0; i < terms.size(); ++i)
  {
    p += terms[i] * monomial;
  }
  return p;
}

/*
  Prints out full polynomial.
*/
template<unsigned int N>
void Polynomial<N>::printP() const
{
  if (terms[0].totalOrder() == 0 || abs(terms[0].c) != 1)
    std::cout << terms[0].c;
  else if (terms[0].c == -1)
    std::cout << "-";
  if (terms[0].c < -1.0e-14 || terms[0].c > 1.0e-14)
    printM(terms[0]);

  for (unsigned int i = 1; i < terms.size(); ++i)
  {
    if (terms[i].c > 0)
      std::cout << " + ";
    else
      std::cout << " - ";
    if (terms[i].totalOrder() == 0 || abs(terms[i].c) != 1)
      std::cout << abs(terms[i].c);
    printM(terms[i]);
  }
  std::cout << std::endl;
}

template<unsigned int N>
bool Polynomial<N>::isZero() const
{
  if (leadingTerm().c > -1.0e-14 && leadingTerm().c < 1.0e-14)
  {
    return true;
  }
  return false;
}

/*
  Sorts terms using bubble sort.  How it's sorted depends
  on the chosen polynomial ordering.
*/
template<unsigned int N>
void Polynomial<N>::sortTerms()
{
  bool sorted = false;
  while (!sorted)
  {
    int swaps = 0;
    for (unsigned int i = 0; i < terms.size() - 1; ++i)
    {
      // Swap elements if they are not in the correct order
      bool needSwap = false;
      std::array<int, N> diff{};
      for (int j = 0; j < N; ++j)
        diff[j] = terms[i][j] - terms[i + 1][j];
      if (ordering == "lex")
      {
        for (int i = 0; i < N; ++i)
        {
          if (diff[i] < 0)
          {
            needSwap = true;
            break;
          }
          else if (diff[i] > 0)
            break;
        }
      }
      else if (ordering == "grlex")
      {
        if (terms[i + 1].totalOrder() > terms[i].totalOrder())
          needSwap = true;
        else if (terms[i + 1].totalOrder() == terms[i].totalOrder())
        {
          for (int i = 0; i < N; ++i)
          {
            if (diff[i] < 0)
            {
              needSwap = true;
              break;
            }
            else if (diff[i] > 0)
              break;
          }
        }
      }
      else if (ordering == "grevlex")
      {
        if (terms[i + 1].totalOrder() > terms[i].totalOrder())
          needSwap = true;
        else if (terms[i + 1].totalOrder() == terms[i].totalOrder())
        {
          for (int i = N - 1; i >= 0; --i)
          {
            if (diff[i] > 0)
            {
              needSwap = true;
              break;
            }
            else if (diff[i] < 0)
              break;
          }
        }
      }

      if (needSwap)
      {
        Term<N> tmp = terms[i];
        terms[i] = terms[i + 1];
        terms[i + 1] = tmp;
        swaps++;
      }
    }
    if (swaps == 0)
      sorted = true;
  }
}

/*
  Adds like terms together and removes terms that are 0
  if the polynomial has n > 1 terms.
*/
template<unsigned int N>
void Polynomial<N>::reduceTerms()
{
  // Adds like terms
  for (unsigned int i = 0; i < terms.size() - 1; ++i)
  {
    if (terms[i + 1].degree() == terms[i].degree())
    {
      terms[i + 1].c += terms[i].c;
      terms.erase(terms.begin() + i);
      --i;
    }
  }
  // Removes terms if they are 0
  for (unsigned int i = 0; i < terms.size(); ++i)
  {
    // if (terms[i].c > -1.0e-14 && terms[i].c < 1.0e-14 && terms.size() > 1)
    if (terms[i].c == 0 && terms.size() > 1)
    {
      terms.erase(terms.begin() + i);
      --i;
    }
  }
}

/*
  Helper function that prints out a single monomial.
*/
template<unsigned int N>
void Polynomial<N>::printM(const Term<N>& t) const
{
  for (int i = 0; i < N; ++i)
  {
    char var = 'z' - N + i + 1;
    if (t(i) > 0)
      std::cout << var;
    if (t(i) > 1)
      std::cout << "^" << t(i);
  }
}

template class Polynomial<3>;
template class Polynomial<26>;