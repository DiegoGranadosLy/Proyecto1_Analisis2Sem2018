/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include <vector>
#include <type_traits>
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <regex>
#include <limits>
#include <algorithm>    // std::sort
#include <vector>
#include <type_traits>
#include "AuxiliaryFunctions.h"
#include <Deflation.hpp>
#include <boost/math/tools/rational.hpp>
using namespace std;
using std::string;
using std::exception;
using std::cout;
using std::abs;
using std::pair;
using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;


namespace anpi {
  
  /// Enumerator makes explicit polishing roots
  enum PolishEnum {
    DoNotPolish,
    PolishRoots
  };

  namespace bmt=boost::math::tools; // for polynomial


    template<class T>
    class Muller {
    private:
        T _xi2 = T(0);
        T _xi1 = T(0);
        T _xi = T(0);
        bool flag = false;
        bool flag1 = false;
        vector <T> solutions = vector<T>();
        const double tolerate2 = 1e-50;
    public:
        Muller(T pFirstX, T pSecondX, T pThirdX) {
          _xi2 = T(pFirstX);
          _xi1 = T(pSecondX);
          _xi = T(pThirdX);


        }

        /*
            Permite llamar a la función dependiendo del tipo de dato que le ingresa, en
            este caso se hace así pues complejos y tipos reales donde ambos tienen funciones
            y operadores diferentes, entonces debe haber un método que permita la implementación de ambos sin copiar tanto código.
            La forma que encontré es usando si el tipo de dato es float y de ahí se obtiene el tipo,
            y se compara con las funciones sobrecargadas donde serán diferenciadas por si es true_type
            o false_type. De esta forma, queda funcionando para cualquier tipo.
            */
        vector <T> findRoots(const polynomial <T> &pPolynomy, bool polish) {
          polynomial <T> polyAux = pPolynomy;
          polynomial <T> residuo;
          while (polyAux.size() > 2 && !flag) {
            flag1 = false;
            T tmp = getRoots(polyAux, typename std::is_floating_point<T>::type());
            if (polish) {
              cout << "POLISH" << endl;
              pulido(tmp, typename std::is_floating_point<T>::type());
              tmp = getRoots(polyAux, typename std::is_floating_point<T>::type());
              //flag1 = true;
            }
            //cout << "La raíz es: "<<tmp<<endl;
            polyAux = deflatefunction(polyAux, tmp, residuo, typename std::is_floating_point<T>::type());
            solutions.push_back(tmp);
            cout << "El polinomio restante es: " << polyAux << endl;
          }
          if (!flag) {
            if (polyAux.size() != 1) solutions.push_back(T(-polyAux[0] / polyAux[1]));
            if (polyAux.size() == 1 && (pPolynomy.size() - 1 != solutions.size())) solutions.push_back(T(-polyAux[0]));
          }

          return solutions;

        }

        void pulido(T solution, std::true_type) {
          if (_xi <= solution) {
            _xi = T(solution);
          } else if (_xi1 <= solution && _xi > solution) {
            _xi1 = T(solution);
          } else if (_xi1 > solution && _xi2 <= solution) {
            _xi2 = T(solution);
          }
        }

        void pulido(T solution, std::false_type) {
          if (_xi.real() < solution.real()) {
            _xi = T(solution.real());
          } else if (_xi1.real() < solution.real() && _xi.real() > solution.real()) {
            _xi1 = T(solution.real());
          } else if (_xi1.real() > solution.real() && _xi2.real() < solution.real()) {
            _xi2 = T(solution.real());
          }
        }

        polynomial <T> deflatefunction(const polynomial <T> &poly, const T &root, polynomial <T> &res, std::true_type) {
          return deflate(poly, root, res);
        }

        polynomial <T>
        deflatefunction(const polynomial <T> &poly, const T &root, polynomial <T> &res, std::false_type) {
          return (abs(root.imag()) == 0) ? deflate(poly, root, res) : deflate2(poly, root, res);
        }

        /*
        Encuentra las raíces de los polinomios si el tipo de dato usado es complejo
        */
        T getRoots(const polynomial <T> pPolynomy, std::false_type) {
          Auxilary<T> *evaluation = new Auxilary<T>();
          T fxi2;
          T fxi1;
          T fxi;
          T q;
          T hSecondX;
          T dFirstX;
          T dSecond;
          T arga;
          T argb;
          T argc;
          T root;
          T _xii2 = _xi2;
          T _xii1 = _xi1;
          T _xii = _xi;
          for (int ite = 0; ite < 100000; ++ite) {
            fxi2 = evaluation->evaluateFunction(pPolynomy, _xii2);
            fxi1 = evaluation->evaluateFunction(pPolynomy, _xii1);
            fxi = evaluation->evaluateFunction(pPolynomy, _xii);
            q = (_xii - _xii1) / (_xii1 - _xii2);
            arga = q * fxi - q * (T(1) + q) * fxi1 + q * q * fxi2;
            argb = (T(2) * q + T(1)) * fxi - (T(1) + q) * (T(1) + q) * fxi1 + q * q * fxi2;
            argc = (T(1) + q) * fxi;
            if ((_xii1 - _xii2) != T(0)) {
              T roottmp = sqrt(argb * argb - T(4) * argc * arga);
              if (argb.real() < 0) root = _xii + (_xii - _xii1) * (T(2) * argc) / (argb - roottmp);
              else root = _xii - (_xii - _xii1) * (T(2) * argc) / (argb + roottmp);

              if ((pow(argb.real(), 2) - 4 * arga.real() * argc.real()) > 0) {
                if (_xii.real() < root.real()) {
                  _xii = T(root.real());
                } else if (_xii1.real() < root.real() && _xii.real() > root.real()) {
                  _xii1 = T(root.real());
                } else if (_xii1.real() > root.real() && _xii2.real() < root.real()) {
                  _xii1 = T(root.real());
                }
              } else {
                _xii2 = T(_xii1.real());
                _xii1 = T(_xii.real());
                _xii = T(root.real());
              }
              T tmpf = evaluation->evaluateFunction(pPolynomy, root);
              if (abs(tmpf.real()) < tolerate2) break;
            } else {
              if ((pow(argb.real(), 2) - 4 * arga.real() * argc.real()) > 0) {
                if (_xi.real() < root.real()) {
                  _xi = T(root.real());
                } else if (_xi1.real() < root.real() && _xi.real() > root.real()) {
                  _xi1 = T(root.real());
                } else if (_xi1.real() > root.real() && _xi2.real() < root.real()) {
                  _xi1 = T(root.real());
                }
              } else {
                _xi2 = T(_xi1.real());
                _xi1 = T(_xi.real());
                _xi = T(root.real());
              }
              break;
            }
          }
          if ((abs(root.imag()) != 0) && !flag1) {
            solutions.push_back(T(root.real(), -root.imag()));
            flag1 = !flag1;
          }
          cout << "Esta raíz es: " << root << endl;
          return root;
        }

        /*
        Encuentra las raíces del polinomio si es dato usado es real.
        */
        T getRoots(const polynomial <T> pPolynomy, std::true_type) {
          Auxilary<T> *evaluation = new Auxilary<T>();
          T fxi2;
          T fxi1;
          T fxi;
          T q;
          T hSecondX;
          T dFirstX;
          T dSecond;
          T arga;
          T argb;
          T argc;
          T root;
          T _xii2 = _xi2;
          T _xii1 = _xi1;
          T _xii = _xi;
          for (int ite = 0; ite < 100; ++ite) {
            fxi2 = evaluation->evaluateFunction(pPolynomy, _xii2);
            fxi1 = evaluation->evaluateFunction(pPolynomy, _xii1);
            fxi = evaluation->evaluateFunction(pPolynomy, _xii);
            q = (_xii - _xii1) / (_xii1 - _xii2);
            arga = q * fxi - q * (T(1) + q) * fxi1 + q * q * fxi2;
            argb = (T(2) * q + T(1)) * fxi - (T(1) + q) * (T(1) + q) * fxi1 + q * q * fxi2;
            argc = (T(1) + q) * fxi;
            if ((_xii1 - _xii2) != 0) {
              T roottmp = sqrt(argb * argb - T(4) * argc * arga);
              if (argb < 0) root = _xii + (_xii - _xii1) * (T(2) * argc) / (argb - roottmp);
              else root = _xii - (_xii - _xii1) * (T(2) * argc) / (argb + roottmp);
              cout << "EL DISCRIMINATE ES: " << (argb * argb - 4 * arga * argc) << endl;
              cout << "LAS Xs SON: " << _xii2 << " " << _xii1 << " " << _xii << endl;
              if ((argb * argb - 4 * arga * argc) < 0) {
                if (!flag) cout << "Este polinomio tiene raices complejas, use otro tipo de dato." << endl;
                flag = !flag;
                _xi2 = T(_xi1);
                _xi1 = T(_xi);
                _xi = T(root);
                break;

              } else {
                if (_xii <= root) {
                  _xii = T(root);
                } else if (_xii1 <= root && _xii > root) {
                  _xii1 = T(root);
                } else if (_xii1 > root && _xii2 <= root) {
                  _xii2 = T(root);
                }
              }
              T tmpf = evaluation->evaluateFunction(pPolynomy, root);
              if (abs(tmpf) < tolerate2) {
                if (_xi <= root) {
                  _xi = T(root);
                } else if (_xi1 <= root && _xi > root) {
                  _xi1 = T(root);
                } else if (_xi1 > root && _xi2 <= root) {
                  _xi2 = T(root);
                }
                break;
              }
            } else {
              if ((pow(argb, 2) - 4 * arga * argc) > 0) {
                if (_xi < root) {
                  _xi = T(root);
                } else if (_xi1 < root && _xi > root) {
                  _xi1 = T(root);
                } else if (_xi1 > root && _xi2 < root) {
                  _xi1 = T(root);
                }
              } else {
                _xi2 = T(_xi1);
                _xi1 = T(_xi);
                _xi = T(root);
              }
              return root;
            }
          }
          return root;
        }
    };

  /**
   * Compute the roots of the given polynomial using the Muller method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
  template<class T,class U>
  void muller(const bmt::polynomial<T>& poly,
              std::vector<U>& roots,
              const PolishEnum polish=DoNotPolish,
              const U start=U()) {
    
    static_assert(std::is_floating_point<T>::value ||
                  boost::is_complex<T>::value,
                  "T must be floating point or complex");
    static_assert(std::is_floating_point<U>::value ||
                  boost::is_complex<U>::value,
                  "U must be floating point or complex");

    throw Exception("Not implemented yet!");
  }
}


#endif
