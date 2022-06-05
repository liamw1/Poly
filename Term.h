#pragma once
#include <array>

/*
  A single polynomial term in R[x_1, x_2,..., x_N].
  Stores the coefficient and the powers on x_1, x_2,..., x_N.
  Terms can be multiplied with "*" or divided with "/"
*/
template<unsigned int N>
class Term
{
public:
  double c = 0;

  Term<N>();

  Term<N>(double iniC, const std::array<unsigned int, N>& iniAlpha);

  unsigned int operator()(unsigned int index) const;

  unsigned int& operator[](unsigned int index);

  Term<N> operator*(const Term<N>& t) const;

  Term<N> operator/(const Term<N>& t) const;

  std::array<unsigned int, N> degree() const;

  unsigned int totalOrder() const;

  bool divides(const Term<N>& t) const;

private:
  std::array<unsigned int, N> alpha{};
};