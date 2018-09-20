/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_DEFLATION_HPP
#define ANPI_DEFLATION_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {
  namespace bmt=boost::math::tools; // bt as alias for boost::math::tools
  
  /*
 * Divide 2 polonomios entre sí
 * @param poly_dividendo es el polinomio dividendo.
 * @param poly_divisor es el polinomio divisor.
 * @param poly_residuo es el polinomio donde se recibe el residuo de la division
 * @return poly_cociente es el polinomio cociente que resulta de la division
 */
template <typename T>
polynomial<T> dividirPolinomios(const polynomial<T>& poly_dividendo, const polynomial<T>& poly_divisor, polynomial<T>& poly_residuo) 
{

                //Grado de los polinomios de entrada
                int grado_dividendo(poly_dividendo.size() - 1);
                int grado_divisor(poly_divisor.size() - 1);

                //Verifica que el polinomio divisor no sea cero
                while (grado_divisor >= 0 && poly_divisor[grado_divisor] == T(0)) {
                    --grado_divisor;
                } if (grado_divisor < 0) throw("Error: dividirPolinomios con divisor igual cero."); 

                //División de polinomios
                poly_residuo = poly_dividendo;
                vector<T> q(grado_dividendo + 1, 0);   //Vector cociente
                for (int k = grado_dividendo - grado_divisor; k >= 0; --k) {
                    q[k] = poly_residuo[k + grado_divisor] / poly_divisor[grado_divisor];   //Genera el cociente
                    for (int j = grado_divisor + k - 1; j >= k; --j) {
                        poly_residuo[j] -= q[k] * poly_divisor[j - k];  //Genera el residuo
                    }
                }
                for (int j = grado_divisor; j <= grado_dividendo; ++j)
                    poly_residuo[j] = T(0); //Completa el residuo con ceros

                polynomial<T>* poly_cociente = new polynomial<T>(q.begin(), q.end());
                return *poly_cociente;
}

template <typename T>
polynomial<T> dividirPolinomiosComplex(const polynomial<T>& poly_dividendo, const polynomial<T>& poly_divisor, polynomial<T>& poly_residuo) 
{

                //Grado de los polinomios de entrada
                int grado_dividendo(poly_dividendo.size() - 1);
                int grado_divisor(poly_divisor.size() - 1);

                //Verifica que el polinomio divisor no sea cero
                while (grado_divisor >= 0 && poly_divisor[grado_divisor].real() == 0) {
                    --grado_divisor;
                } if (grado_divisor < 0) throw("Error: dividirPolinomios con divisor igual cero."); 

                poly_residuo = poly_dividendo;
                vector<T> q(grado_dividendo + 1, 0);

                for (int k = grado_dividendo - grado_divisor; k >= 0; --k) {
                    q[k] = poly_residuo[k + grado_divisor] / poly_divisor[grado_divisor];
                    for (int j = grado_divisor + k - 1; j >= k; --j) {
                        poly_residuo[j] -= q[k] * poly_divisor[j - k];
                    }
                }
                for (int j = grado_divisor; j <= grado_dividendo; ++j)
                    poly_residuo[j] = T(0);

                polynomial<T>* poly_cociente = new polynomial<T>(q.begin(), q.end());
                return *poly_cociente;
}
  /**
   * Deflate polynomial
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
  template<class T>
  bmt::polynomial<T> deflate(const bmt::polynomial<T>& poly, const T& root, T& residuo, T& tolerance=anpi::epsilon<T>())
  {
    int const grado_poly(poly.size() - 1);     
    //Crear polinomio a partir de raiz
    //Crear x - raiz
    int grado_poly_root(1);
    vector<T> vpoly_root(grado_poly_root + 1);
    vpoly_root[1] = T(1);
    vpoly_root[0] = T(-root);
    polynomial<T> poly_root(vpoly_root.begin(), vpoly_root.end());
    return dividirPolinomios(poly, poly_root, residuo);    
  }

  /**
   * Deflate polynomial with a second order polynomial.
   *
   * The second order polynomial equals x^2 -2 Re(root)x + |root|^2.
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
  template<class T>
  bmt::polynomial<T> deflate2(const bt::polynomial<T>& poly, const std::complex<T>& root, bt::polynomial<T>& residuo, T& tolerance=anpi::epsilon<T>())
  {
  //Declara variables locales
  int grado_poly(poly.size() - 1);
  T root_real(real(root));
  T root_imag(imag(root));           
  //Crear polinomio a partir de raiz compleja
  //Crear x^2 - 2ax + (a^2 + b^2)
  int grado_poly_root_c(2);
  vector<T> poly_root_c(grado_poly_root_c + 1);
  poly_root_c[2] = T(1);
  poly_root_c[1] = T(-2)*root_real;
  poly_root_c[0] = T(root_real*root_real + root_imag*root_imag);
  polynomial<T> poly_root(poly_root_c.begin(), poly_root_c.end());
  return dividirPolinomios(poly, poly_root, residuo);
  } 
}



#endif
