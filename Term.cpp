#include "Term.h"
#include <iostream>

template<unsigned int N>
Term<N>::Term() {}

template<unsigned int N>
Term<N>::Term(double iniC, const std::array<unsigned int, N>& iniAlpha)
{
  c = iniC;
  alpha = iniAlpha;
}

template<unsigned int N>
unsigned int Term<N>::operator()(unsigned int index) const
{
  return alpha[index];
}

template<unsigned int N>
unsigned int& Term<N>::operator[](unsigned int index)
{
  return alpha[index];
}

template<unsigned int N>
Term<N> Term<N>::operator*(const Term<N>& t) const
{
  std::array<unsigned int, N> gamma;
  for (int i = 0; i < N; ++i)
    gamma[i] = alpha[i] + t(i);
  return Term(c * t.c, gamma);
}

template<unsigned int N>
Term<N> Term<N>::operator/(const Term<N>& t) const
{
  std::array<unsigned int, N> gamma;
  for (int i = 0; i < N; ++i)
    gamma[i] = alpha[i] - t(i);
  return Term(c / t.c, gamma);
}

template<unsigned int N>
std::array<unsigned int, N> Term<N>::degree() const
{
  return alpha;
}

template<unsigned int N>
unsigned int Term<N>::totalOrder() const
{
  unsigned int totalOrder = 0;
  for (int i = 0; i < N; ++i)
    totalOrder += alpha[i];
  return totalOrder;
}

template<unsigned int N>
bool Term<N>::divides(const Term<N>& t) const
{
  for (int i = 0; i < N; ++i)
  {
    if (alpha[i] > t(i))
      return false;
  }
  return true;
}

template class Term<3>;
template class Term<26>;