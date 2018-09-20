/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

/**
 * El código fuente presentado a continuación fue tomado del proyecto
 * simbody, el cual cuenta con derechos de uso como se detalla más abajo,
 * y se le realizaron algunas ligeras modificaciones y adaptaciones para
 * el proyecto. Recurso disponible en https://github.com/simbody/simbody/tree/master/SimTKcommon/Polynomial/src
 */
/*
*******************************************************************************
*
*
*                       Copyright (c) 2002
*                       Henrik Vestermark
*                       Denmark
*
*                       All Rights Reserved
*
*   Permission to use, copy, distribute, and sell this software and its
*   documentation for any purpose is hereby granted without fee, provided:
*   THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
*   EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY WARRANTY
*   OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL
*   Henrik Vestermark or Future Team Aps, BE LIABLE FOR ANY SPECIAL,
*   INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES
*   WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER OR NOT
*   ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF LIABILITY,
*   ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
*   SOFTWARE.
*
*******************************************************************************
*
*
* Module name     :   cpoly.cpp
* Module ID Nbr   :
* Description     :   cpoly.cpp -- Jenkins-Traub real polynomial root finder.
*                     Translation of TOMS493 from FORTRAN to C. This
*                     implementation of Jenkins-Traub partially adapts
*                     the original code to a C environment by restruction
*                     many of the 'goto' controls to better fit a block
*                     structured form. It also eliminates the global memory
*                     allocation in favor of local, dynamic memory management.
*
*                     The calling conventions are slightly modified to return
*                     the number of roots found as the function value.
*
*                     INPUT:
*                     opr - vector of real coefficients in order of
*                          decreasing powers.
*                     opi - vector of imaginary coefficients in order of
*                          decreasing powers.
*                     degree - integer degree of polynomial
*
*                     OUTPUT:
*                     zeror,zeroi - output vectors of the
*                          real and imaginary parts of the zeros.
*                            to be consistent with rpoly.cpp the zeros is inthe index
*                            [0..max_degree-1]
*
*                     RETURN:
*                     returnval:   -1 if leading coefficient is zero, otherwise
*                          number of roots found.
* --------------------------------------------------------------------------
* Change Record   :
*
* Version  Author/Date     Description of changes
* -------  -----------     ----------------------
* 01.01    HVE/20021101      Initial release
* 01.02    PE/20070808       Converted to a class, templatized
* 01.03    MAS/20130410      Minor changes to match RPoly changes
*
* End of Change Record
* --------------------------------------------------------------------------
*/

#ifndef ANPI_JENKINS_TRAUB_HPP
#define ANPI_JENKINS_TRAUB_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include "math.h"
#include <float.h>

namespace anpi
{
    namespace bmt=boost::math::tools; // for polynomial
    using namespace std;
    using namespace boost::math::tools;

// Defining every single function
	template <typename T>
	static void noshft(const int l1, const int nn, T h[], T p[], const double eta, T t);

	template <typename T>
	static void fxshft(const int l2, T *z, int *conv, const int nn, T s, T p[], T qp[], T *pv, T *t, T h[], T sh[], T qh[], const double are, const double mre, const double eta, const double infin);

	template <typename T>
	static void vrshft(const int l3, T *z, int *conv, T s, T p[], T qp[], T *pv, const int nn, const double are, const double mre, const double eta, const double infin, T h[], T qh[], T *t);

	template <typename T>
	static void calct(int *bol, const int nn, T s, T h[], T qh[], const double are, T *pv, T *t);

	template <typename T>
	static void nexth(const int bol, const int nn, T qh[], T qp[], T *t, T h[]);

	template <typename T>
	static void polyev(const int nn, const T s, const T p[], T q[], T *pv);

	template <typename T>
	static double errev(const int nn, const T complejo[], const double ms, const double mp, const double are, const double mre);

	template <typename T>
	static void cauchy(const int nn, T pt[], double *fn_val);

	template <typename T>
	static double scale(const int nn, const T pt[], const double eta, const double infin, const double smalno, const double base);

	template <typename T>
	static void cdivid(const T a, const T b, T *c);

	template <typename T>
	static double cmod(const T complejo);

	template <class T, class U>
	typename bmt::polynomial<std::complex<typename anpi::detail::inner_type<U>::type>> castCoeffToResult(bmt::polynomial<T> poly);

	static void mcon(double *eta, double *infiny, double *smalno, double *base);


    // MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE PROGRAM.
    // THE USER MAY EITHER SET THEM DIRECTLY OR USE THE STATEMENTS BELOW TO
    // COMPUTE THEM. THE MEANING OF THE FOUR CONSTANTS ARE -
    // ETA       THE MAXIMUM RELATIVE REPRESENTATION ERROR WHICH CAN BE DESCRIBED
    //           AS THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    //           1.0_dp + ETA &gt; 1.0.
    // INFINY    THE LARGEST FLOATING-POINT NUMBER
    // SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
    // BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED
    //
    static void mcon(double *eta, double *infiny, double *smalno, double *base)
    {
        *base = numeric_limits<double>::radix;
        *eta = numeric_limits<double>::epsilon();
        *infiny = numeric_limits<double>::max();
        *smalno = numeric_limits<double>::min();
    }


    // RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
    // THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
    // INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
    // BASE.
    // PT - MODULUS OF COEFFICIENTS OF P
    // ETA, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.
    template <typename T>
    static double scale(const int nn, const T pt[], const double eta, const double infin, const double smalno, const double base)
    {
        int i, l;
        double hi, lo, max, min, x, sc;
        double fn_val;

        // Find largest and smallest moduli of coefficients
        hi = sqrt(infin);
        lo = smalno / eta;
        max = 0;
        min = infin;

        for (i = 0; i <= nn; i++)
        {
            x = pt[i].real();
            if (x > max)
                max = x;
            if (x != 0 && x < min)
                min = x;
        }

        // Scale only if there are very large or very small components
        fn_val = 1;
        if (min >= lo && max <= hi)
            return fn_val;
        x = lo / min;
        if (x <= 1)
            sc = 1 / (sqrt(max) * sqrt(min));
        else
        {
            sc = x;
            if (infin / sc > max)
                sc = 1;
        }
        l = (int)(log(sc) / log(base) + 0.5);
        fn_val = pow(base, l);
        return fn_val;
    }

    // COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
    //
    template <typename T>
    static void cdivid(const T a, const T b, T *c)
    {
        double r, d, t, infin;

        if (b.real() == 0 && b.imag() == 0)
        {
            // Division by zero, c = infinity
            mcon(&t, &infin, &t, &t);
            c->real(infin);
            c->imag(infin);
            return;
        }

        if (fabs(b.real()) < fabs(b.imag()))
        {
            r = b.real() / b.imag();
            d = b.imag() + r * b.real();
            c->real((a.real() * r + a.imag()) / d);
            c->imag((a.imag() * r - a.real()) / d);
            return;
        }

        r = b.imag() / b.real();
        d = b.real() + r * b.imag();
        c->real((a.real() + a.imag() * r) / d);
        c->imag((a.imag() - a.real() * r) / d);
    }

    // CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
    // POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
    template<typename T>
    static void cauchy( const int nn, T pt[], double *fn_val )
    {
        int i, n;
        double x, xm, f, dx, df;

        pt[nn].real(-pt[nn].real());

        // Compute upper estimate bound
        n = nn;
        x = exp(log(-pt[nn].real()) - log(pt[0].real())) / n;
        if( pt[n-1].real() != 0 )
        {
            // Newton step at the origin is better, use it
            xm = -pt[nn].real() / pt[ n - 1 ].real();
            if(xm < x) x = xm;
        }

        // Chop the interval (0,x) until f < 0
        while(1)
        {
            xm = x * 0.1;
            f = pt[0].real();
            for(i = 1; i <= nn; i++)
                f = f * xm + pt[ i ].real();
            if(f <= 0)
                break;
            x = xm;
        }
        dx = x;

        // Do Newton iteration until x converges to two decimal places
        while(fabs( dx / x ) > 0.005)
        {
            pt[ 0 ].imag(pt[ 0 ].real());
            for(i = 1; i <= nn; i++)
                pt[i].imag(pt[i-1].imag() * x + pt[i].real());
            f = pt[nn].imag();
            df = pt[0].imag();
            for(i = 1; i < n; i++)
                df = df * x + pt[ i ].imag();
            dx = f / df;
            x -= dx;
        }

        *fn_val = x;
    }

    // COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
    // POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
    template<typename T>
    static void noshft(const int l1, const int nn, T h[], T p[], const double eta, T t)
    {
        int i, j, jj, n, nm1;
        double xni, t1, t2;
        n = nn;
        nm1 = n - 1;
        for(i = 0; i < n; i++)
        {
            xni = nn - i;
            h[i].real(xni * p[i].real() / n);
            h[i].imag(xni * p[i].imag() / n);
        }
        for(jj = 1; jj <= l1; jj++)
        {
            if(cmod<T>( h[n-1]) > eta * 10 * cmod<T>(p[n-1]))
            {
                cdivid<T>(-p[nn], h[n-1], &t);
                for(i = 0; i < nm1; i++)
                {
                    j = nn - i - 1;
                    t1 = h[j-1].real();
                    t2 = h[j-1].imag();
                    h[j].real(t.real() * t1 - t.imag() * t2 + p[j].real());
                    h[j].imag(t.real() * t2 + t.imag() * t1 + p[j].imag());
                }
                h[0].real(p[0].real());
                h[0].imag(p[0].imag());
            }
            else
            {
                // If the constant term is essentially zero, shift H coefficients
                for(i = 0; i < nm1; i++)
                {
                    j = nn - i - 1;
                    h[j].real(h[j-1].real());
                    h[j].imag(h[j-1].imag());
                }
                h[0].real(0);
                h[0].imag(0);
            }
        }
    }

    // COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
    // INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
    // APPROXIMATE ZERO IF SUCCESSFUL.
    // L2 - LIMIT OF FIXED SHIFT STEPS
    // ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
    // CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
    //
    template<typename T>
    static void fxshft(const int l2, T *z, int *conv, const int nn, T s, T p [], T qp [], T *pv, T *t, T h[], T sh[],T qh[], const double are, const double mre, const double eta,const double infin)
    {
        int i, j, n, test, pasd, bol;
        T ot, svs;
        n = nn;
        polyev<T>(nn, s, p, qp, pv);
        test = 1;
        pasd = 0;

        // Calculate first T = -P(S)/H(S)
        //int *bol, const int nn, T s, T h[], T qh[], const double are, T pv, T* t
        calct<T>(&bol,nn,s,h,qh,are,pv,t);

        // Main loop for second stage
        for(j = 1; j <= l2; j++)
        {
            ot.real(t->real());
            ot.imag(t->imag());

            // Compute the next H Polynomial and new t
            //const int bol, const int nn, T qh [], T qp[], T t, T h[]
            nexth<T>(bol,nn,qh,qp,t,h);
            calct<T>( &bol,nn,s,h,qh,are,pv,t );
            z->real(s.real() + t->real());
            z->imag(s.imag() + t->imag());

            // Test for convergence unless stage 3 has failed once or this
            // is the last H Polynomial
            T temp;
            temp.real(t->real() - ot.real());
            temp.imag(t->imag() - ot.imag());
            if(!(bol || !test || j == 12)){
                if(cmod<T>(temp) < 0.5 * cmod<T>((*z)))
                {
                    if(pasd)
                    {
                        // The weak convergence test has been passwed twice, start the third stage
                        // Iteration, after saving the current H polynomial and shift
                        for(i = 0; i < n; i++)
                        {
                            sh[i].real(h[i].real());
                            sh[i].imag(h[i].imag());
                        }
                        svs.real(s.real());
                        svs.imag(s.imag());
                        // const int l3, T *z, int *conv, T s, T p[], T qp [], T* pv, const int nn, const double are, const double mre, const double eta,const double infin, T h[], T qh[], T *t
                        //const int l2, T *z, int *conv, const int nn, T s, T p [], T qp [], T pv, T *t, T h[], T sh[],T qh[], const double are, const double mre, const double eta,const double infin
                        vrshft<T>(10,z,conv,s,p,qp,pv,nn,are,mre,eta,infin,h,qh,t);
                        if(*conv) return;

                        //The iteration failed to converge. Turn off testing and restore h,s,pv and T
                        test = 0;
                        for(i = 0; i < n; i++)
                        {
                            h[i].real(sh[i].real());
                            h[i].imag(sh[i].imag());
                        }
                        s.real(svs.real());
                        s.imag(svs.imag());
                        polyev<T>(nn, s, p, qp, pv);
                        calct<T>(&bol, nn, s, h, qh,are,pv,t);
                        continue;
                    }
                    pasd = 1;
                }
                else
                    pasd = 0;
            }
        }

        // Attempt an iteration with final H polynomial from second stage
        vrshft<T>(10, z, conv,s,p,qp,pv,nn,are,mre,eta,infin,h,qh,t);
    }

    // CARRIES OUT THE THIRD STAGE ITERATION.
    // L3 - LIMIT OF STEPS IN STAGE 3.
    // ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
    //           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
    // CONV    -  .TRUE. IF ITERATION CONVERGES
    //
    template<typename T>
    static void vrshft( const int l3, T *z, int *conv, T s, T p[], T qp [], T* pv, const int nn, const double are, const double mre, const double eta,const double infin, T h[], T qh[], T *t)
    {
        int b, bol, i, j;
        double mp, ms, omp, relstp, r1, r2, tp;

        *conv = 0;
        b = 0;
        s.real(z->real());
        s.imag(z->imag());

        // Main loop for stage three
        for(i = 1; i <= l3; i++)
        {
            // Evaluate P at S and test for convergence
            polyev<T>(nn, s, p, qp, pv);
            mp = cmod<T>(*pv);
            ms = cmod<T>(s);
            if(mp <= 20 * errev<T>(nn,qp,ms,mp,are,mre))
            {
                // Polynomial value is smaller in value than a bound onthe error
                // in evaluationg P, terminate the ietartion
                *conv = 1;
                z->real(s.real());
                z->imag(s.imag());
                return;
            }
            if(i != 1)
            {
                if(!(b || mp < omp || relstp >= 0.05))
                {
                    // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed
                    // shift steps into the cluster to force one zero to dominate
                    tp = relstp;
                    b = 1;
                    if(relstp < eta) tp = eta;
                    r1 = sqrt(tp);
                    r2 = s.real() * (1 + r1) - s.imag() * r1;
                    s.imag(s.real() * r1 + s.imag() * (1 + r1));
                    s.real(r2);
                    polyev<T>(nn, s, p, qp, pv);
                    for(j = 1; j <= 5; j++)
                    {
                        calct<T>(&bol,nn,s,h,qh,are,pv,t);
                        nexth<T>(bol,nn,qh,qp,t,h);
                    }
                    omp = infin;
                    goto _20;
                }

                // Exit if polynomial value increase significantly
                if(mp *0.1 > omp) return;
            }

            omp = mp;

            // Calculate next iterate
            _20:  calct<T>(&bol,nn,s,h,qh,are,pv,t);
            nexth<T>(bol,nn,qh,qp,t,h);
            calct<T>(&bol,nn,s,h,qh,are,pv,t);
            if(!bol)
            {
                relstp = cmod<T>(*t) / cmod<T>(s);
                s.real(s.real()+t->real());
                s.imag(s.imag()+ t->imag());
            }
        }
    }

    // COMPUTES  T = -P(S)/H(S).
    // BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
    template<typename T>
    static void calct( int *bol, const int nn, T s, T h[], T qh[], const double are, T *pv, T* t)
    {
        int n;
        T hv;

        n = nn;

        // evaluate h(s)
        polyev<T>(n - 1, s, h, qh, &hv);
        *bol = cmod<T>(hv) <= are * 10 * cmod<T>(h[n-1]) ? 1 : 0;
        if(!*bol)
        {
            cdivid<T>(-*pv, hv, t);
            return;
        }

        t->real(0);
        t->imag(0);
    }

    // CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
    // BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
    //
    template<typename T>
    static void nexth( const int bol, const int nn, T qh [], T qp[], T *t, T h[] )
    {
        int j, n;
        double t1, t2;

        n = nn;
        if(!bol)
        {
            for(j = 1; j < n; j++)
            {
                t1 = qh[j-1].real();
                t2 = qh[j-1].imag();
                h[j].real(t->real() * t1 - t->imag() * t2 + qp[j].real());
                h[j].imag(t->real() * t2 + t->imag() * t1 + qp[j].imag());
            }
            h[0].real(qp[0].real());
            h[0].imag(qp[0].imag());
            return;
        }

        // If h[s] is zero replace H with qh
        for(j = 1; j < n; j++)
        {
            h[j].real(qh[j-1].real());
            h[j].imag(qh[j-1].imag());
        }
        h[0].real(0);
        h[0].imag(0);
    }

    // EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
    // PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
    //
    template<typename T>
    static void polyev( const int nn, const T s, const T p[], T q[], T *pv )
    {
        int i;
        double t;

        q[0].real(p[0].real());
        q[0].imag(p[0].imag());
        pv->real(q[0].real());
        pv->imag(q[0].imag());

        for( i = 1; i <= nn; i++ )
        {
            t = (pv->real()) * s.real() - (pv->imag()) * s.imag() + p[i].real();
            pv->imag((pv->real()) * s.imag() + (pv->imag()) * s.real() + p[i].imag());
            pv->real(t);
            q[i].real(pv->real());
            q[i].imag(pv->imag());
        }
    }

    // BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
    // QR,QI - THE PARTIAL SUMS
    // MS    -MODULUS OF THE POINT
    // MP    -MODULUS OF POLYNOMIAL VALUE
    // ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION
    //
    template<typename T>
    static double errev( const int nn, const T complejo[], const double ms, const double mp, const double are, const double mre )
    {
        int i;
        double e;

        e = cmod<T>(complejo[0]) * mre / (are + mre);
        for(i = 0; i <= nn; i++)
            e = e * ms + cmod<T>( complejo[i]);

        return e * (are + mre) - mp * mre;
    }

    // MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
    template <typename T>
    static double cmod(const T complejo)
    {
      T a;

      a.real(fabs(complejo.real()));
      a.imag(fabs(complejo.imag()));
      if (a.real() < a.imag())
        return a.imag() * sqrt(1.0 + pow((a.real() / a.imag()), 2.0));

      if (a.real() > a.imag())
        return a.real() * sqrt(1.0 + pow((a.imag() / a.real()), 2.0));

      return (double)(a.real() * sqrt(2.0));
    }

    template<typename T>
    int cpoly(polynomial <T> op, int degree, vector <T> &zero)
    {
        static double are, mre, eta, infin;
        static int nn;
        static T sComp, tComp, pvComp, zComp;
        static T *pComp, *hComp, *qpComp, *qhComp, *shComp;
        int cnt1, cnt2, idnn2, i, conv;
        double xx, yy, cosr, sinr, smalno, base, xxx, bnd;

        mcon(&eta, &infin, &smalno, &base);
        are = eta;
        mre = 2.0 * sqrt(2.0) * eta;
        xx = 0.70710678;
        yy = -xx;
        cosr = -0.060756474;
        sinr = -0.99756405;
        nn = degree;

        // Algorithm fails if the leading coefficient is zero, or degree = 0
        if (degree == 0 && op[0].real() == 0 && op[0].imag() == 0)
        {
            return -1;
        }

        // Remove the zeros at the origin if any
        while(op[nn].real() == 0 && op[nn].imag() == 0)
        {
            idnn2 = degree - nn;
            zero[idnn2].real(0);
            zero[idnn2].imag(0);
            nn--;
        }

        if (nn == 0)
            return degree;

        // Allocate arrays
        pComp = new T [degree+1];
        hComp = new T [degree+1];
        qpComp= new T [degree+1];
        qhComp = new T [degree+1];
        shComp= new T [degree+1];

        // Make a copy of the coefficients
        for(i = 0; i <= nn; i++)
        {
            pComp[i].real(op[i].real());
            pComp[i].imag(op[i].imag());
            shComp[i].real(cmod<T>( pComp[i]));
        }

        // Scale the polynomial
        bnd = scale<T>(nn,shComp,eta, infin, smalno, base );
        if(bnd != 1)
            for(i = 0; i <= nn; i++)
            {
                pComp[i].real(pComp[i].real()*bnd) ;
                pComp[i].imag(pComp[i].imag() *bnd);
            }

        search:
        if(nn <= 1)
        {
            cdivid<T>(-pComp[1], pComp[0], &zero[degree-1]);
            delete[] pComp;
            delete[] hComp;
            delete[] qpComp;
            delete[] qhComp;
            delete[] shComp;
            return degree;
        }

        for(i = 0; i<= nn; i++)
            shComp[i].real(cmod<T>(pComp[i]));

        cauchy<T>(nn, shComp, &bnd);
        // Outer loop to control 2 Major passes with different sequences of shifts
        for(cnt1 = 1; cnt1 <= 2; cnt1++)
        {
            // First stage  calculation , no shift
            // const int l1, const int nn, T h[], T p[], const double eta, T t
            noshft<T>(5, nn, hComp, pComp, eta , tComp);
            // Inner loop to select a shift
            for(cnt2 = 1; cnt2 <= 9; cnt2++)
            {
                // Shift is chosen with modulus bnd and amplitude rotated by 94 degree from the previous shif
                xxx = cosr * xx - sinr * yy;
                yy = sinr * xx + cosr * yy;
                xx = xxx;
                sComp.real(bnd * xx);
                sComp.imag(bnd * yy);
                // Second stage calculation, fixed shift
                fxshft<T>( 10 * cnt2, &zComp, &conv ,nn ,sComp,pComp, qpComp,&pvComp,&tComp,hComp, shComp,qhComp, are,mre,eta,infin);
                if(conv)
                {
                    // The second stage jumps directly to the third stage ieration
                    // If successful the zero is stored and the polynomial deflated
                    idnn2 = degree - nn;
                    zero[idnn2].real(zComp.real());
                    zero[idnn2].imag(zComp.imag());
                    nn--;
                    for(i = 0; i <= nn; i++)
                    {
                        pComp[i].real(qpComp[i].real());
                        pComp[i].imag(qpComp[i].imag());
                    }
                    goto search;
                }
                // If the iteration is unsuccessful another shift is chosen
            }
            // if 9 shifts fail, the outer loop is repeated with another sequence of shifts
        }
        // The zerofinder has failed on two major passes
        // return empty handed with the number of roots found (less than the original degree)
        degree -= nn;
        // Deallocate arrays
        delete [] pComp;
        delete [] hComp;
        delete [] qpComp;
        delete [] qhComp;
        delete [] shComp;

        return degree;
    }

  /**
   * Compute the roots of the given polynomial using the Jenkins-Traub method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
  template<class T,class U>
  void jenkinsTraub(const bmt::polynomial<T>& poly,
                    std::vector<U>& roots) {
    
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
