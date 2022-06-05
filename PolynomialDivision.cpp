#include <iostream>
#include <chrono>
#include "Term.h"
#include "Polynomial.h"

#define MAX(a,b) ((a > b) ? a : b)

struct Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<float> duration;
  const char* timerName;

  Timer(const char* name) : timerName(name)
  {
    duration = std::chrono::high_resolution_clock::duration::zero();
  }

  void timeStart()
  {
    start = std::chrono::high_resolution_clock::now();
  }

  void timeEnd()
  {
    end = std::chrono::high_resolution_clock::now();
    duration += end - start;
  }

  ~Timer()
  {
    std::cout << timerName << " took " << duration.count() * 1000.0f << "ms" << std::endl;
  }
};

/*
  Prints out the quotients and remainder at some point in the division algorithm.
  Set step = -1 to print final results.
*/
template<unsigned int N>
void printStep(int step, const Polynomial<N>& p, const std::vector<Polynomial<N>>& Q, const Polynomial<N>& r)
{
  if (step != -1)
  {
    std::cout << "Step " << step << ":" << std::endl;
    std::cout << "p: ";
    p.printP();
  }
  else std::cout << "Final Results:" << std::endl;
  for (unsigned int i = 0; i < Q.size(); ++i)
  {
    std::cout << "q" << i + 1 << ": ";
    Q[i].printP();
  }
  std::cout << "r: ";
  r.printP();
  std::cout << std::endl;
  if (step == -1)
  {
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }
}

/*
  Returns remainder of f upon division by a set of Polynomials F.
*/
template<unsigned int N>
Polynomial<N> MultivariatePolynomialDivision(Polynomial<N> f, const std::vector<Polynomial<N>>& F, const bool printEveryStep)
{
  unsigned int s = F.size();
  const char* ordering = f.getOrdering();

  // Create list of quotients and remainder
  std::vector<Polynomial<N>> Q;
  for (unsigned int i = 0; i < s; ++i)
    Q.emplace_back(Polynomial<N>(ordering));
  Polynomial<N> r = Polynomial<N>(ordering);

  Polynomial<N>& p = f;

  // Division algorithm
  int step = 1;
  while (p.leadingTerm().c != 0)
  {
    unsigned int i = 0;
    bool divisionOccurred = false;
    while (i < s && !divisionOccurred)
    {
      // Checks if the leading term of f_i divides the leading term of p
      if (F[i].leadingTerm().divides(p.leadingTerm()))
      {
        Q[i] += p.leadingTerm() / F[i].leadingTerm();
        p -= F[i] * (p.leadingTerm() / F[i].leadingTerm());
        divisionOccurred = true;

        if (printEveryStep)
        {
          printStep(step, p, Q, r);
          step++;
        }
      }
      else
      {
        i++;
      }
    }
    if (!divisionOccurred)
    {
      r += p.leadingTerm();
      p -= p.leadingTerm();

      if (printEveryStep)
      {
        printStep(step, p, Q, r);
        step++;
      }
    }
  }

  // Print final results
  if (printEveryStep)
  {
    printStep(-1, p, Q, r);
  }

  return r;
}

/*
  Returns the least common multiple of the leading monomials of f and g.
*/
template<unsigned int N>
Term<N> LeadingMonomialLCM(const Polynomial<N>& f, const Polynomial<N>& g)
{
  std::array<unsigned int, N> alpha = f.leadingTerm().degree();
  std::array<unsigned int, N> beta = g.leadingTerm().degree();
  std::array<unsigned int, N> gamma{};
  for (int i = 0; i < N; ++i)
    gamma[i] = MAX(alpha[i], beta[i]);
  return Term<N>(1.0, gamma);
}

/*
  Computes the S-Polynomial of f and g.
*/
template<unsigned int N>
Polynomial<N> S_Polynomial(const Polynomial<N>& f, const Polynomial<N>& g)
{
  return f * (LeadingMonomialLCM(f, g) / f.leadingTerm()) - g * (LeadingMonomialLCM(f, g) / g.leadingTerm());
}

/*
  Converts a set of Polynomials F into a Grobner basis using Buchberger's Algorithm.
*/
template<unsigned int N>
void ConvertToGrobnerBasis(std::vector<Polynomial<N>>& F)
{
LOOP:
  std::vector<Polynomial<N>> G = F;
  for (unsigned int i = 0; i < G.size(); ++i)
  {
    for (unsigned int j = 0; j < G.size(); ++j)
    {
      if (i != j)
      {
        Polynomial<N> r = MultivariatePolynomialDivision(S_Polynomial(G[i], G[j]), G, false);
        if (!r.isZero())
        {
          F.emplace_back(r);
          goto LOOP;
        }
      }
    }
  }
}