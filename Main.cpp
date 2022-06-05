#include "PolynomialDivision.cpp"


int main()
{
  // Problem 3a
  const char* ordering = "lex";
  const int n = 3;

  std::vector<Polynomial<n>> G = { Polynomial<n>("x + y + z - 3", ordering), Polynomial<n>("x^2 + y^2 + z^2 - 5", ordering), Polynomial<n>("x^3 + y^3 + z^3 - 7", ordering) };
  ConvertToGrobnerBasis(G);

  for (unsigned int i = 0; i < G.size(); ++i)
    G[i].printP();
  std::cout << std::endl;

  Polynomial<n> f = Polynomial<n>("x^4 + y^4 + z^4 - 9", ordering);
  MultivariatePolynomialDivision(f, G, false).printP();
  std::cout << std::endl;


  // Problem 3b
  Polynomial<n> g1 = Polynomial<n>("x^5 + y^5 + z^5", ordering);
  MultivariatePolynomialDivision(g1, G, false).printP();

  Polynomial<n> g2 = Polynomial<n>("x^6 + y^6 + z^6", ordering);
  MultivariatePolynomialDivision(g2, G, false).printP();
}