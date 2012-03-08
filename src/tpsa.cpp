/*
 * Copyright(C) 2008 by Lingyun Yang
 * lingyun(.dot.]yang@gmail.com
 * http://www.lingyunyang.com
 *
 * Please get permission from Lingyun Yang before you redistribute this file.
 *
 */

//! \brief Automatic Differentiation Test
//! \file tpsa.cpp
//! \version $Id: tpsa.cpp,v 1.4 2009-06-08 10:48:43 frs Exp $
//! \author Lingyun Yang, http://www.lingyunyang.com/

#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "tpsa.h"

static unsigned int comb_num(unsigned int n, unsigned int r);
static unsigned int gcd(unsigned int a, unsigned int b);
static void init_order_index(unsigned int nv, unsigned int nd);
static void init_prod_index(unsigned int nv, unsigned int nd);
static bool within_limit(TNVND nv, TNVND nd);
static TNVND* choose(TNVND n, TNVND r, TNVND* p, TNVND nv, TNVND nd);
static void init_base(TNVND nv, TNVND nd);
static void print_vec(unsigned int, std::ostream&);

//! Search for the element index of the product of two elements.
// static unsigned int prod_index(TNVND* pb, TNVND* pm, TNVND order);

//! Maximum length of AD vector. A system with nv variables and nd orders will
//! have length C(nv+nd, nd).
static const unsigned int MAX_N_BASE =
    std::numeric_limits<unsigned int>::max()/2;

//! global nv, nd
static TNVND gnv, gnd;

//! C(nv+nd, nd)
static unsigned int FULL_VEC_LEN;

//! Automatic Differentiation
//! starting index of elements with certain orders.
static unsigned int* order_index;


static unsigned int** H;

//! base vector.
//! e.g. (000) to (004) if 4th order 3 variables, C(nv+nd,nd)
//! vectors in total.
static TNVND *base;
static unsigned int** prdidx; // index of product
static TNVND nvec; // number of vec assigned.
static TNVND* vec; // vec initialized.
// static unsigned int tblsize; // not used

// memory pool
static double **advecpool;
static std::vector<double*> advec;
static std::vector<unsigned int> adveclen;


using namespace std;

//! Greatest Common Divisor, Euclidean algorithm.
unsigned int gcd(unsigned int a, unsigned int b)
{
    unsigned int t;
    while (b != 0u) {
        t = b;
        b = a % b;
        a = t;
    }
    return a;
}

//! \brief Combinatorial Number (Binomial Coefficient, choose r from n)
//
// The result should fit into a unsigned int(e.g. on 32bit system 4294967296).
//
// C(15,10) = 3003
// C(17,10) = 19448
// C(18,10) = 43758
// int16    = 65536
// C(19,10) = 92378
// C(20,10) = 184756
// C(30,10) = 30045015
// int32    = 4294967296
unsigned int comb_num(unsigned int n, unsigned int r)
{
    const unsigned int MAX_N_R = 100;

    if (n == 0 || r == 0) return 1;

    // not allowed, regardless of mathematical meaning.
    if (n < r) return 0;

    // obviously too large fitting in a unsigned int
    if (r > MAX_N_R && (n-r) > MAX_N_R) return 0;

    unsigned int k = (r > n-r ? n-r : r);

    std::vector<int> numerator(k), denominator(k);

    for (unsigned int i = 0; i < k; ++i) {
        numerator[i] = n - i;
        denominator[i] = k - i;
    }

    int c = 1;
    for (size_t i = 0; i < k; ++i) {
        // deal with denominator[i]
        for (size_t j = 0; j < k; ++j) {
            if ( (c=gcd(denominator[i], numerator[j])) > 1 ) {
                if (numerator[j] % c != 0 || denominator[i] % c != 0) {
                 #ifdef DEBUG
                    std::cerr << "ERROR: gcd error, " << k << " = gcd("
                              << denominator[i] << ","
                              << numerator[j] << ")" << std::endl;
                 #endif

                    return 0;
                }
                numerator[j] /= c;
                denominator[i] /= c;
            }
        }
    }

    // check if the result is an integer, i.e. the denominator vector is
    // (1,1,1, ... ,1)
    for (size_t i = 0; i < k; ++i) {
        if (denominator[i] != 1) {
         #ifdef DEBUG
            std::cerr << "ERROR: denominator is not dividable by numerator"
                      << std::endl;
         #endif

            return 0;
        }
    }

    size_t N = 1;

    for (size_t i = 0; i < k; ++i) {
        if (MAX_N_BASE*1.0/numerator[i] < N) {
            return 0;
        }
        N *= numerator[i];
    }
    return N;
}

//! Initialize the index of orders up to nv+1.
//! init the index of each order in monimial vec
void init_order_index(unsigned int nv, unsigned int nd)
{
    order_index = new unsigned int[nd+2];
    size_t i = 0, N = 0;
    order_index[0] = 0;

    while (i < nd+1) {
        N += comb_num(nv+i-1, i);
        order_index[++i] = N;
    }
    order_index[nd+1] = FULL_VEC_LEN;
}

//! Initialize hash table
void init_prod_index(unsigned int nv, unsigned int nd)
{
    //unsigned int m[nv+1][nd+2];
    unsigned int** m;
    m = new unsigned int* [nv+1];
    for (size_t i = 0; i < nv+1; ++i) {
        m[i] = new unsigned int[nd+2];
    }

    TNVND* vtmp = new TNVND[nv+1];
    // hash table
    H = new unsigned int*[nv+1];

    vtmp[0] = 0;
    for (TNVND i = 0; i < nv+1; ++i) {
        H[i] = new unsigned int[nd+2];

        m[i][0] = H[i][0] = 0;
        m[i][1] = H[i][1] = 1;
        for (TNVND j = 2; j < nd+2; ++j) {
            m[i][j] = H[i][j] = m[i][j-1]*(i+j-2)/(j-1);
        }
        for (TNVND j = 1; j < nd+2; ++j) {
            m[i][j] += m[i][j-1];
            H[i][j] += H[i][j-1];
        }
    }

 #ifdef PRT_H
    int width_h = 3;
    if (order_index[nd] > 999) width_h = 5;
    else if (order_index[nd] > 99) width_h = 4;

    std::cout << std::setw(width_h) << ' ';
    for (TNVND i = 0; i < nv + 2; ++i) {
        std::cout << std::setw(width_h) << (unsigned int)i;
    }
    std::cout << std::endl;
    std::cout << "   ------------------" << std::endl;

    for (TNVND i = 0; i < nv+1; ++i) {
        std::cout << std::setw(4) << (unsigned int) i;
        for (TNVND j = 0; j < nd+2; ++j)
            std::cout << std::setw(4) << m[i][j];
        std::cout << std::endl;
    }
    std::cout << std::endl;
 #endif

    unsigned NS = 0;
    prdidx = new unsigned int*[FULL_VEC_LEN];
    prdidx[0] = NULL;
    TNVND ord = 1;
    TNVND *pb=base;

    // allocating memory.
    for (size_t i = 1; i < order_index[nd]; ++i) {
        if (order_index[ord+1] <= i) ++ord;
        size_t M = order_index[nd-ord+1];
        prdidx[i] = new unsigned int[M];
        NS += M;
        for (size_t j = 1; j < M; ++j) {
            prdidx[i][j] = 0;
            size_t idx = 0;
            for(TNVND k = 0; k < nv; ++k) {
                idx += m[nv-k][*(pb+i*nv+k) + *(pb+j*nv+k)];
            }
            prdidx[i][j] = idx;
        }
    }
    // tblsize = NS; not used
    //std::cout << "Mem: " << NS << std::endl;
}

//! \brief Reserve space for n number of TPS vectors.
//!
/** \param n number of TPS vectors
 * \see ad_init
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_reserve(const TVEC* n)
#else
void ad_reserve(const unsigned int* n)
#endif
{
    const unsigned int N = *n;

    if (N <= 0) return;

    advecpool = new double*[N];
    for (size_t i = 0; i < N; ++i) {
        advecpool[i] = new double[FULL_VEC_LEN];
        advec.push_back(NULL);
        if (adveclen.size() <= i) adveclen.push_back(0);
        else adveclen[i] = 0;
    }

 #ifdef DEBUG_ALL
    std::cout << "Initialize: " << advec.size() << " vectors" << std::endl;
 #endif
}



//! \brief Check if C(nd+nv, nv) is in the maximum number
//! allowed. i.e. C(nv+nd, nv) < MAX_N_BASE.
bool within_limit(TNVND nv, TNVND nd)
{
    if (nv == 0 || nd == 0) return true;

    unsigned int N = comb_num(nv+nd, nv);

    // too large, can not fit in an unsigned int.
    if (N == 0) return false;

    if ( N >= MAX_N_BASE ) return false;
    else return true;
}

//! \brief Get all the combinations, choose r from n.
//! Do nothing if n==0 or r==0. other numbers are not checked.
TNVND* choose(TNVND n, TNVND r, TNVND* p, TNVND nv, TNVND /* nd */)
{
    std::vector<TNVND> c(r);
    TNVND i, j, k;

    // size_t n = 0;
    // Boundary condition
    if (n == 0 || r == 0) return 0;

    for (i = 0; i < r; ++i) {
        c[i] = i;
    }
    c[r-1] = r-2;
    while(true) {
        j = r - 1;
        // for unsigned types, this is tricky.
        // two more constraint j < r, c[j] < n.
        // j >= 0 is always true
        while(j < r && c[j] >= n - r +j && c[j] < n) j--;
        // if (j < 0) break;
        if (j >= r ) break;

        ++c[j];
        for ( k = j + 1; k < r; ++k) c[k] = c[k-1] + 1;

        /*
          for (TNVND iv = 0; iv < r; ++iv)
          std::cout << ' ' << static_cast<unsigned int>(c[iv]-iv);
          // std::cout << std::endl;
          std::cout << "     ";
        */

        // CoeffMatrix.push_back(c);
        for (size_t iv = 0; iv < r; ++iv)
            ++(*(p+c[iv]-iv));
        p += nv;


        /*
          TNVND* h = p - nv;
          for (size_t iv = 0; iv < nv; ++iv) {
          std::cout << ' ' << static_cast<unsigned int>(*h);
          ++h;
          }
          std::cout << std::endl;
        */

        // ++n;
    }

    return p;
}

//
void init_base(TNVND nv, TNVND nd)
{
    // int nv = variable; // max num of variables
    // unsigned int nm = 0;
    // tps_base.push_back(base);
    base = new TNVND[nv*FULL_VEC_LEN];
    for (size_t i = 0; i < nv*FULL_VEC_LEN; ++i) base[i] = 0;

    TNVND* pb = base;

    if ( !within_limit(nv, nd) ) {
        return;
    }

    // std::cout << "Address: " << (unsigned int)(base) << std::endl;
    for (size_t i = 0; i < nv; ++i) {
        *pb = 0;
        ++pb;
    }

    // std::cout << "Address: " << (unsigned int)(pb) << std::endl;
    for (size_t d = 1; d <= nd; ++d) {
        pb = choose(nv+d-1, d, pb, nv, nd);
        // std::cout << "Address: " << (unsigned int)(pb) << std::endl;

    }
    // std::cout << "Total: " << pb << std::endl;
    pb = base;
    for (size_t i = 0; i < FULL_VEC_LEN; ++i) {
        TNVND x = 0;
        for (TNVND k = 0; k < nv; ++k) {
            x += *(pb+nv-1-k);
            *(pb+nv-1-k) = x;
        }
        pb += nv;
    }
}

/**\brief Initialize the TPSA library.
 *
 * \param nvar number of variables
 * \param nord highest order of TPS.
 *
 * This must be called before any use of TPS vector, and followed by ad_reserve
 */
#ifdef MSVC_DLL
extern "C" _declspec(dllexport) void _stdcall ad_init(const TNVND* nvar, const TNVND* nord)
#else
void ad_init(const TNVND* nvar, const TNVND* nord)
#endif
{
    TNVND nv = *nvar;
    TNVND nd = *nord;

 #ifdef DEBUG
    std::cerr << "Initialize nv= " << *nvar << " nd= " << *nord << std::endl;
 #endif
    advecpool = NULL;
    gnv = nv; // assign to global nv
    gnd = nd; // assign to global nd

    FULL_VEC_LEN = comb_num(gnv+gnd, gnd);
    init_order_index(gnv, gnd);

    // overflow !!
    if (FULL_VEC_LEN == 0) {
        std::cerr << "Overflow!  (" << nv << "," << nd << ")" << std::endl;
        std::exit(1);
    }

    // this reset the base vector
    vec = new TNVND[nv];
    for (TNVND i = 0; i < nv; ++i) vec[i] = 0;

    nvec = 0;

    init_base(gnv, gnd);

    init_prod_index(gnv, gnd);

    // random number generator
    srand(1U);
}

/** \brief Reset the base vectors.
 * To avoid confusions when starting a new set of TPSA calculations, call this before a new ad_init().
 * with nv from ad_nvar(nv).
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_resetvars(const TNVND* nvar)
#else
void ad_resetvars(const TNVND* nvar)
#endif
{
    if (!vec) return;
    if (*nvar > gnv) {
        //std::cerr << "ERROR: Reset more base vectors than earlier initialization. "
        //          << std::endl;
        //std::cerr << "   global nv = " << gnv << "  input= " << *nvar << std::endl;
        for (TNVND i = 0; i < gnv; ++i) vec[i] = 0;
    } else {
        for (TNVND i = 0; i < *nvar; ++i) vec[i] = 0;
    }
    nvec = 0;
}

/** \brief Allocate space for one TPS vector.
 * \param i return a new integer representing the allocated TPS vector.
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_alloc(TVEC* i)
#else
void ad_alloc(unsigned int* i)
#endif
{
    // didn't check the upper limit of i
    *i = advec.size();
    for (size_t j = 0; j < advec.size(); ++j) {
        if (!advec[j]) {
            *i = j;
            break;
        }
    }

    if (*i >= advec.size()) {
        std::cerr << "Run out of vectors" << std::endl;
        exit(-1);
        //return;
    }


    //size_t n = 0;
    //ad_count(&n);
    //std::cerr << "Found space for AD vector: " << *i << " " << n << std::endl;
    unsigned int iv = *i;

    // initialize as a constant 0, length 1
    advec[iv] = advecpool[iv];
    for (size_t j = 0; j < FULL_VEC_LEN; ++j) {
        advec[iv][j] = 0;
    }
    adveclen[iv] = 1;
    advec[iv][0] = 0;

 #ifdef DEBUG_ALL
    std::cout << "Allocate: " << iv << " len: " << adveclen[iv] << std::endl;
    ad_print(&iv);
 #endif
}

/** \brief Free the memory allocated by ad_alloc
 * \param i TPS vector
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_free(const TVEC* i)
#else
void ad_free(const unsigned int* i)
#endif
{
    advec[*i] = NULL;
    adveclen[*i] = 0;

 #ifdef DEBUG_ALL
    std::cout << "AD free " << *i << std::endl;
 #endif
}

/** \brief Copy TPS vector
 * \param isrc source
 * \param idst destination
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_copy(const TVEC* isrc, const TVEC* idst)
#else
void ad_copy(const unsigned int* isrc, const unsigned int* idst)
#endif
{
    unsigned int i = *isrc;
    unsigned int j = *idst;
    if (i == j) return;
    // did not check if j is allocated.
    for (size_t k = 0; k < FULL_VEC_LEN; ++k)
        advec[j][k] = advec[i][k];
    adveclen[j] = adveclen[i];
 #ifdef DEBUG_ALL
    std::cout << "Copy from " << i << " to " << j << std::endl;
 #endif
    //ad_print(idst);
}

/** \brief Count how many vectors has been allocated
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_count(TVEC* n)
#else
void ad_count(TVEC* n)
#endif
{
    *n = 0;
    for(size_t i = 0; i < advec.size(); ++i)
        if (advec[i]) ++(*n);
}

/** \brief Set a TPS vector as a constant
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_const(const TVEC* iv, const double* r)
#else
void ad_const(const TVEC* iv, const double* r)
#endif
{
    unsigned int ii = *iv;
    for(size_t i = 0; i < adveclen[ii]; ++i)
        advec[ii][i] = 0;
    advec[ii][0] = *r;
 #ifdef DEBUG_ALL
    std::cout << "Set const: iv " << *iv << " = " << *r << std::endl;
 #endif
    adveclen[ii] = 1;
}

/** \brief Set the small TPS coefficients as 0 if they are less than eps.
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_clean(const TVEC* iv, const double* eps)
#else
void ad_clean(const TVEC* iv, const double* eps)
#endif
{
    unsigned int N = 0;
    for(size_t i = 0; i < adveclen[*iv]; ++i) {
        if (abs(advec[*iv][0]) < abs(*eps)) advec[*iv][i] = 0;
        else N = i;
    }
    if (adveclen[*iv] > N+1) {
        adveclen[*iv] = N+1;
    }
}

/** \brief Fill the TPS coefficients as a random number in [0,1]
 *
 * The length(highest nonzero coefficients) is not changed.
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm)
#else
void ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm)
#endif
{
    for(size_t i = 0; i < adveclen[*iv]; ++i) {
        if (std::rand()*1.0/RAND_MAX > *ratio)
            advec[*iv][i] = 0;
        else advec[*iv][i] = std::rand()*(*xm)/RAND_MAX;
        //std::cout << advec[*iv][i] << " ";
    }
    //std::cout << std::endl;
}

/**
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_pok(const TVEC* ivec, int* c, size_t* n, double* x)
#else
void ad_pok(const unsigned int* ivec, int* c, size_t* n, double* x)
#endif
{
    const unsigned int ii = *ivec;
    size_t N = (*n > gnv) ? gnv : *n;
    double r = *x;
 #ifdef DEBUG_ALL
    std::cout << "AD Pok ivec= " << *ivec << " n= " << *n << " x= " << *x << std::endl;
 #endif

    unsigned int *cef = new unsigned int[gnv];
    unsigned int *bv = new unsigned int[gnv];
    unsigned int d = 0;
    for (size_t i = 0; i < N; ++i) {
        //std::cout << " " << c[i];
        cef[i] = c[i];
        d += c[i];
    }
    // didn't check if
    for (size_t i = N; i < gnv; ++i) cef[i] = 0;

    //std::cout  << std::endl;
    //for (int i = 0; i < FULL_VEC_LEN; ++i) {
    //    advec[ii][i] = 0;
    //}

    if (d > gnd) return;

    size_t k = 0;

    //std::cout << "order: " << (unsigned int)d << std::endl;
    for (size_t i = 0; i < gnv; ++i){
        bv[i] = d;
        d -= cef[i];
        //std::cout << (unsigned int)i << " : "
        //          << (unsigned int) t << ' '
        //          << (unsigned int)bv[i] << ' '
        //          << (unsigned int)d << std::endl;
        k += H[gnv-i][bv[i]];
    }
    //std::cout << std::endl;
    //std::cout << k << std::endl;
    advec[ii][k] = r;

    if (k+1>adveclen[ii])
        adveclen[ii] =  k+1;
    delete []cef;
    delete []bv;
}

// ith element of ivec, save pattern in c, value in x
// c has length gnv
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x)
#else
void ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x)
#endif
{
    const TVEC ii = *ivec;

    for (size_t i = 0; i < gnv; ++i) c[i] = 0;
    if (*idx > adveclen[ii] || *idx < 1) {
        * x = 0;
        return;
    }

    //
    double* v = advec[ii];
    *x = v[*idx-1];
    TNVND* p = base + gnv*(*idx-1);
    for (size_t j = 0; j < gnv-1; ++j) {
        c[j] = (unsigned int) (*p-*(p+1));
        ++p;
    }
    c[gnv-1] = (unsigned int) *p;
}


#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_pek(const TVEC* ivec, int* c, size_t* n, double* x)
#else
void ad_pek(const unsigned int* ivec, int* c, size_t* n, double* x)
#endif
{
    const unsigned int ii = *ivec;
    size_t N = (*n > gnv) ? gnv : *n;
    //double r = *x;

 #ifdef DEBUG_ALL
    std::cout << "AD Pek ivec= " << *ivec << " n= " << *n << " x= " << *x << std::endl;
 #endif

    unsigned int *cef = new unsigned int[gnv];
    unsigned int *bv = new unsigned int[gnv];
    unsigned int d = 0;
    for (size_t i = 0; i < N; ++i) {
        //std::cout << " " << c[i];
        cef[i] = c[i];
        d += c[i];
    }
    // didn't check if
    for (size_t i = N; i < gnv; ++i) cef[i] = 0;

    //std::cout  << std::endl;
    //for (int i = 0; i < FULL_VEC_LEN; ++i) {
    //    advec[ii][i] = 0;
    //}

    if (d > gnd) return;

    size_t k = 0;

    //std::cout << "order: " << (unsigned int)d << std::endl;
    for (size_t i = 0; i < gnv; ++i){
        bv[i] = d;
        d -= cef[i];
        //std::cout << (unsigned int)i << " : "
        //          << (unsigned int) t << ' '
        //          << (unsigned int)bv[i] << ' '
        //          << (unsigned int)d << std::endl;
        k += H[gnv-i][bv[i]];
    }
    //std::cout << std::endl;
    //std::cout << k << std::endl;

    if (k > adveclen[ii]) *x = 0;
    else *x = advec[ii][k];

    delete []bv;
    delete []cef;
}

// set the vector ii be the base vector iv
// ii, iv are 0-started index.
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_var(const TVEC* ivec, const double* x, unsigned int* ibvec)
#else
void ad_var(const unsigned int* ivec, const double* x, unsigned int* ibvec)
#endif
{
    const unsigned int iv = *ivec;
    unsigned int ibv = *ibvec;
    double x0 = *x;
    //unsigned int nbvmax;
    // TNVND i = iv - 1;
    double *v = advec[iv];
 #ifdef DEBUG_ALL
    std::cout << "Enter ad_var " << *ivec << "  " << *x << "  " << *ibvec << std::endl;
 #endif

    for (size_t k = 0; k < FULL_VEC_LEN; ++k) v[k] = 0;
    v[0] = x0;
    //v[iv] = 1;

    // check iv and is not set before
    if (ibv < gnv && vec[ibv] == 0) {
        ++vec[ibv];
        ++nvec;
        adveclen[iv] = ibv + 2; // length is const+first order, 1+(iv+1)
        v[ibv+1] = 1;
    } else if (ibv >= gnv ) {
        std::cerr << "Out of boundary, init as an ordinary variable"
                  << std::endl;
        adveclen[iv] = 1;
    } else if (vec[ibv] > 0)  {
      #ifdef DEBUG
        // this will update the vec[iv]
        if (vec[ibv] == 1) {
            std::cerr << "WARNING: Base vector " << (unsigned int) ibv
                      << " has been initialized [" << gnv << "]: (";
            for (TNVND j = 0; j < gnv; ++j) {
                if (vec[j]) std::cerr << ' ' << (unsigned int)vec[j];
                else std::cerr << " N";
            }
            std::cerr << ")" << std::endl;
        }
        ++vec[ibv];
    #endif
        adveclen[iv] = ibv + 2;
        v[ibv+1] = 1;
    } else {
        std::cerr << "What else ?" << std::endl;
    }

 #ifdef DEBUG_ALL
    std::cout << "Exit ad_var" << std::endl;
    ad_print(&iv);
 #endif
    //std::exit(-1);
    //for (size_t k = 1; k < iv; ++k) v[k] = 0;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_abs(const TVEC* iv, double* r)
#else
void ad_abs(const TVEC* iv, double* r)
#endif
{
    *r = 0;
    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        *r += std::abs(advec[*iv][i]);
        //std::cout << advec[*iv][i] << ' ';
    }
    //std::cout << "abs: " << *r << std::endl;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_truncate(const TVEC* iv, const TNVND* d)
#else
void ad_truncate(const TVEC* iv, const TNVND* d)
#endif
{
    if (*d > gnd) return;

    for(size_t i = order_index[*d]; i < adveclen[*iv]; ++i) {
        advec[*iv][i] = 0;
    }
    adveclen[*iv] = order_index[*d];
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_add(const TVEC* idst, const TVEC* jsrc)
#else
void ad_add(const unsigned int* idst, const unsigned int* jsrc)
#endif
{
    unsigned int i = *idst;
    unsigned int j = *jsrc;

    double *v = advec[i];
    double *rhsv = advec[j];
    //double *resv = advec[k];

    if (adveclen[i] < adveclen[j]) {
        // two blocks for overlap part and non-overlap part.
        for (size_t ii = 0; ii < adveclen[i]; ++ii) v[ii] += rhsv[ii];
        for (size_t ii = adveclen[i]; ii < adveclen[j]; ++ii) v[ii] = rhsv[ii];
        adveclen[i] = adveclen[j];
    }else {
        for (size_t ii = 0; ii < adveclen[j]; ++ii) v[ii] += rhsv[ii];
    }

 #ifdef DEBUG_ALL
    std::cout << "AD add " << *idst << " " << *jsrc << std::endl;
 #endif

    //while(std::abs(advec[*idst][adveclen[*idst]-1]) < std::numeric_limits<double>::min()  && adveclen[*idst] > 2) --adveclen[*idst];

    //return *this;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_add_const(const TVEC* i, double *r)
#else
void ad_add_const(const TVEC* i, double *r)
#endif
{
    advec[*i][0] += *r;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_sub(const TVEC* idst, const TVEC* jsrc)
#else
void ad_sub(const unsigned int* idst, const unsigned int* jsrc)
#endif
{
    unsigned int i = *idst;
    unsigned int j = *jsrc;

    double *v = advec[i];
    double *rhsv = advec[j];

 #ifdef DEBUG_ALL
    std::cout << "AD sub " << *idst << " ["  << adveclen[i]
              << "]  " << *jsrc << "  [" << adveclen[j] << "]"
              << std::endl;
 #endif

    //double *resv = advec[k];
    if (adveclen[i] == 0 || adveclen[j] == 0) {
        std::cerr << "ERROR: AD sub zero length vector" << std::endl;
        //std::exit(-1);
        return;
    }

    if (adveclen[i] < adveclen[j]) {
        // two blocks for overlap part and non-overlap part.
        for (size_t ii = 0; ii < adveclen[i]; ++ii) v[ii] -= rhsv[ii];
        for (size_t ii = adveclen[i]; ii < adveclen[j]; ++ii) v[ii] = -rhsv[ii];
        adveclen[i] = adveclen[j];
    }else {
        for (size_t ii = 0; ii < adveclen[j]; ++ii) v[ii] -= rhsv[ii];
    }

    //while(abs(advec[*idst][adveclen[*idst]-1]) < std::numeric_limits<double>::min()  && adveclen[*idst] > 1) --adveclen[*idst];

    //return *this;
}

//! internal multiplication, dst should be different from lhs and rhs.
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_mult(const TVEC* ilhs, const TVEC* irhs, TVEC* idst)
#else
void ad_mult(const unsigned int* ilhs, const unsigned int* irhs, unsigned int* idst)
#endif
{
    unsigned int lhs = *ilhs;
    unsigned int rhs = *irhs;
    unsigned int dst = *idst;

 #ifdef DEBUG_ALL
    std::cout << "AD mult " << *ilhs << " " << *irhs
              << " " << *idst << std::endl;
    ad_print(ilhs);
    ad_print(irhs);
 #endif

    // can speed up without copy-constructor.
    for (size_t i = 0; i < adveclen[dst]; ++i) advec[dst][i] = 0;

    adveclen[dst] = adveclen[lhs];

    advec[dst][0] = advec[lhs][0] * advec[rhs][0];
    for (size_t i = 1; i < adveclen[rhs]; ++i) {
        advec[dst][i] += advec[lhs][0]*advec[rhs][i];
    }
    for (size_t i = 1; i < adveclen[lhs]; ++i) {
        advec[dst][i] += advec[lhs][i]*advec[rhs][0];
    }

 #ifdef DEBUG_MULT
    int width_v = 1;
    if (nv > static_cast<TNVND>(9)) width_v = 2;
 #endif

    // It is symmetric, but not when lhs rhs have different size, .....
    TNVND ord = 1;
    size_t L = std::max(adveclen[lhs], adveclen[rhs]);
 #ifdef DEBUG_ALL
    std::cerr << "Length: " << L << " nv " << gnv << "  nd " << gnd << std::endl;
 #endif

    for (size_t i = 1; i < std::min(adveclen[lhs], order_index[gnd]); ++i) {
        if (order_index[ord+1] <= i) ++ord;
        unsigned int M = order_index[gnd-ord+1];
        if (M > adveclen[rhs]) M = adveclen[rhs];
        //dst.v[prdidx[i][i]] += lhs.v[i]*rhs.v[i];
        //if (prdidx[i][i] > L) L = prdidx[i][i];
        for (size_t j = 1; j < M; ++j) {
            advec[dst][prdidx[i][j]] += advec[lhs][i]*advec[rhs][j];
            //dst.v[prdidx[i][j]] += lhs.v[j]*rhs.v[i];
            //if (prdidx[i][j] >= L) L = prdidx[i][j] + 1;
            //std::cout << i << ',' << j << "  " << L << std::endl;
        }
        //std::cout << i << ',' << M-1 << "  p=" << prdidx[i][M-1] << " L=" << L;
        if (prdidx[i][M-1] >= L) {
            L = prdidx[i][M-1] + 1;
            //std::cout << ' ' << L << std::endl;
        }
    }

    //while(abs(dst.v[L-1]) < std::numeric_limits<T>::min()) --L;
    if (L > FULL_VEC_LEN) L = FULL_VEC_LEN;
    adveclen[dst] = L;

    while(std::abs(advec[dst][adveclen[dst]-1]) < std::numeric_limits<double>::min()  && adveclen[dst] > 1) --adveclen[dst];

    //std::cerr << "Len: " << adveclen[dst] << std::endl;
    //ad_print(idst);
    //dst.len = FULL_VEC_LEN;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_mult_const(const TVEC* iv, double* c)
#else
void ad_mult_const(const TVEC* iv, double* c)
#endif
{
    double* v = advec[*iv];
    for (size_t i = 0; i < adveclen[*iv]; ++i)
        v[i] *= (*c);
 #ifdef DEBUG_ALL
    std::cout << "AD mult const " << *iv << " const= " << *c << std::endl;
 #endif
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_div_c(const TVEC* iv, const double* c)
#else
void ad_div_c(const TVEC* iv, const double* c)
#endif
{
    if (std::abs(*c) < std::numeric_limits<double>::min()) {
        std::cerr << "ERROR: divide a two small number! " << *c << std::endl;
        std::exit(-1);
    }

    for (size_t i = 0; i < adveclen[*iv]; ++i)
        advec[*iv][i] /= *c;
}

//
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_c_div(const TVEC* iv, const double* c, TVEC* ivret)
#else
void ad_c_div(const TVEC* iv, const double* c, TVEC* ivret)
#endif
{
    TVEC ipn, ip, itmp;
    ad_alloc(&ipn);
    ad_alloc(&ip);
    //ad_alloc(&iret);
    ad_alloc(&itmp);

    TVEC iret = *ivret;

    double *pn = advec[ipn];
    double *p = advec[ip];
    double *ret = advec[iret];
    double *v = advec[*iv];

    ad_copy(iv, &ip);
    ad_copy(iv, &ipn);

    // TODO: didn't check c = 0.0
    double x0 = v[0];
 #ifdef DEBUG_ALL
    std::cout << "x0 " << x0 << std::endl;
 #endif
    for (size_t i = 1; i < adveclen[ip]; ++i) {
        p[i] /= -x0;
        pn[i] = p[i];
    }
    pn[0] = p[0] = 0;

    ret[0] = 1;
    adveclen[iret] = 1;
    //ret += p;
 #ifdef DEBUG_ALL
    std::cout << "-- ret[0]= " << ret[0] << " " << adveclen[iret]
              << " " << p[0] << " " << adveclen[ip] << std::endl;
 #endif

    ad_add(&iret, &ip);
    //ad_print(&ip);
 #ifdef DEBUG_ALL
    std::cout << "-- ret[0]= " << ret[0] << " " << adveclen[iret]
              << " " << p[0] << " " << adveclen[ip]<< std::endl;
 #endif

    for (TNVND nd = 2; nd < gnd+1; ++nd) {
        ad_mult(&ipn, &ip, &itmp);
        ad_copy(&itmp, &ipn);
      for (size_t i = 0; i < adveclen[ipn]; ++i)
          ret[i] += pn[i];
     #ifdef DEBUG_ALL
        std::cout << "ret[0]= " << ret[0] << std::endl;
     #endif
    }
    adveclen[iret] = FULL_VEC_LEN;

    //std::cout << "AD " << iret << std::endl;
    for (size_t i = 0; i < adveclen[iret]; ++i) {
        ret[i] /= x0/(*c);
        //std::cout << ' ' << ret[i];
    }
    //std::cout << std::endl;
 #ifdef DEBUG_ALL
    std::cout << ret[0] << std::endl;
 #endif

    ad_free(&itmp);
    ad_free(&ip);
    ad_free(&ipn);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst)
#else
void ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst)
#endif
{
    TVEC itmp;
    double c = 1.0;
    bool nonzero = false;
    for (size_t i = 0; i < adveclen[*irhs]; ++i) {
        if(std::abs(advec[*irhs][i]) > std::numeric_limits<double>::min()) {
            nonzero = true;
            break;
        }
    }

    if (!nonzero) {
        std::cerr << "ERROR: Divided by zero: " << std::endl;
        ad_print(irhs);
        c = std::sqrt(-1.0);
        std::cerr << c << std::endl;
        std::exit(-1);
    }

    ad_alloc(&itmp);
    ad_c_div(irhs, &c, &itmp);

 #ifdef DEBUG_ALL
    std::cout << "AD div" << std::endl;
    ad_print(&itmp);
    std::cout << "AD mult" << std::endl;
    ad_print(ilhs);
 #endif

    ad_mult(ilhs, &itmp, idst);

    //ad_print(idst);
    ad_free(&itmp);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_sqrt(const TVEC* iv, const TVEC* iret)
#else
void ad_sqrt(const TVEC* iv, const TVEC* iret)
#endif
{
    double x = advec[*iv][0];
    double c;
    TVEC itmp, ip, ipn;
    ad_alloc(&itmp);
    ad_alloc(&ip);
    ad_alloc(&ipn);
    ad_copy(iv, &ip);
    ad_div_c(&ip, &x);
    advec[ip][0] = 0;
    ad_reset(iret);
    advec[*iret][0] = 1;
    adveclen[*iret] = 1;
    // ip is the delta_x
    ad_copy(&ip, &itmp);
    ad_copy(&ip, &ipn);
    c = .5;
 #ifdef DEBUG_ALL
    std::cout << "Coef: c= " << c << std::endl;
 #endif
    for(size_t i = 1; i < gnd+1; ++i) {
        ad_mult_const(&itmp, &c);
        ad_add(iret, &itmp);
        c = c*(1-2.0*i)/2/(i+1.0);
     #ifdef DEBUG_ALL
        std::cout << "   " << i << ":" << c << std::endl;
     #endif
        ad_mult(&ip, &ipn, &itmp);
        ad_copy(&itmp, &ipn);
    }
 #ifdef DEBUG_ALL
    std::cout << std::endl;
 #endif
    x =std::sqrt(x);
    ad_mult_const(iret, &x);
    ad_free(&ipn);
    ad_free(&ip);
    ad_free(&itmp);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_exp(const TVEC* iv, const TVEC* iret)
#else
void ad_exp(const TVEC* iv, const TVEC* iret)
#endif
{
    double x = std::exp(advec[*iv][0]);
    double c;
    TVEC itmp, ip, ipn;

 #ifdef DEBUG_ALL
    std::cout << "AD exp " << *iv << " " << x << "  ret: " << *iret << std::endl;
 #endif

    ad_alloc(&itmp);
    ad_alloc(&ip);
    ad_alloc(&ipn);

    ad_copy(iv, &ip);
    //ad_div_c(&ip, &x);
    advec[ip][0] = 0;

    ad_reset(iret);
    advec[*iret][0] = 1.0;
    adveclen[*iret] = 1;
    // ip is the delta_x

 #ifdef DEBUG_ALL
    std::cout << "reset" << std::endl;
 #endif

    ad_copy(&ip, &itmp);
    ad_copy(&ip, &ipn);
    c = 1.0;
 #ifdef DEBUG_ALL
    std::cout << "Coef: c= " << c << std::endl;
 #endif

    for(size_t i = 1; i < gnd+1; ++i) {
        c = c*1.0*i;
        ad_div_c(&itmp, &c);
        ad_add(iret, &itmp);
     #ifdef DEBUG_ALL
        std::cout << "   " << i << ":" << c << std::endl;
     #endif
        ad_mult(&ip, &ipn, &itmp);
        ad_copy(&itmp, &ipn);
    }
 #ifdef DEBUG_ALL
    std::cout << std::endl;
 #endif
    ad_mult_const(iret, &x);
    ad_free(&ipn);
    ad_free(&ip);
    ad_free(&itmp);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_log(const TVEC* iv, const TVEC* iret)
#else
void ad_log(const TVEC* iv, const TVEC* iret)
#endif
{
    double x = std::log(advec[*iv][0]);
    double c;
    TVEC itmp, ip, ipn;

 #ifdef DEBUG_ALL
    std::cout << "AD exp " << *iv << " " << x << "  ret: " << *iret << std::endl;
 #endif

    ad_alloc(&itmp);
    ad_alloc(&ip);
    ad_alloc(&ipn);

    ad_copy(iv, &ip);
    ad_div_c(&ip, advec[*iv]);
    advec[ip][0] = 0;

    ad_reset(iret);
    advec[*iret][0] = x;
    adveclen[*iret] = 1;
    // ip is the delta_x

 #ifdef DEBUG_ALL
    std::cout << "reset" << std::endl;
 #endif

    ad_copy(&ip, &itmp);
    ad_copy(&ip, &ipn);
    c = 1.0;
 #ifdef DEBUG_ALL
    std::cout << "Coef: c= " << c << std::endl;
 #endif

    for(size_t i = 1; i < gnd+1; ++i) {
        c = (i % 2 == 0 ? -1.0*i: (1.0*i) );
        ad_div_c(&itmp, &c);
        ad_add(iret, &itmp);
     #ifdef DEBUG_ALL
        std::cout << "   " << i << ":" << c << std::endl;
     #endif

        ad_mult(&ip, &ipn, &itmp);
        ad_copy(&itmp, &ipn);
    }
 #ifdef DEBUG_ALL
    std::cout << std::endl;
 #endif
    //ad_mult_const(iret, &x);
    ad_free(&ipn);
    ad_free(&ip);
    ad_free(&itmp);
}

/** \brief Reset to a constant 0
 */
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_reset(const TVEC* iv)
#else
void ad_reset(const TVEC* iv)
#endif
{
 #ifdef DEBUG_ALL
    std::cout << "AD reset " << *iv << std::endl;
 #endif

    for (size_t i = 0; i < adveclen[*iv]; ++i)
        advec[*iv][i] = 0;
    adveclen[*iv] = 0;
}


#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_sin(const TVEC* iv, const TVEC* iret)
#else
void ad_sin(const TVEC* iv, const TVEC* iret)
#endif
{
    //AD<T,nv,nd> ret(v), pnev(v), pnod(v), p(v);
    TVEC ipnev, ipnod, ip;
    ad_alloc(&ipnev);
    ad_alloc(&ipnod);
    ad_alloc(&ip);

    ad_copy(iv, iret);
    ad_copy(iv, &ipnev);
    ad_copy(iv, &ipnod);
    ad_copy(iv, &ip);

    double *v = advec[*iv];
    double *ret = advec[*iret];
    double *pnod = advec[ipnod];
    double *pnev = advec[ipnev];
    double *p = advec[ip];

 #ifdef DEBUG_ALL
    std::cout << ' ' << *iret << ' ' << ret
              << ", " << ipnev << ' ' << pnev
              << ", " << ipnod << ' ' << pnod
              << ", " << ip << ' ' << p
              << std::endl;
 #endif
    // TODO: didn't check c = 0.0
    double s = sin(v[0]);
    double c = cos(v[0]);

    // odd and even
//    size_t pnevlen = adveclen[ipnev], pnodlen=adveclen[ipnod]; // not used

    pnev[0] = pnod[0] = p[0] = 0;

    ret[0] = s;

    for (size_t i = 1; i < adveclen[ip]; ++i) {
        ret[i] *= c;
    }

    //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
    for (TNVND k = 2; k < gnd+1; ++k) {
        //std::cout << (unsigned int) p.v << std::endl << p << std::endl;
        //std::cout << (unsigned int) pnod.v << std::endl << pnod << std::endl;
        //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
        ad_mult(&ip, &ipnod, &ipnev);
        //AD<T,nv,nd>::mult(p, pnod, pnev);
        //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
      for (size_t i = 0; i < adveclen[ipnev]; ++i) {
            pnev[i] /= 1.0*k;
            switch(k % 4) {
            case 0:
                ret[i] += s*pnev[i];
                break;
            case 1:
                ret[i] += c*pnev[i];
                break;
            case 2:
                ret[i] -= s*pnev[i];
                break;
            case 3:
                ret[i] -= c*pnev[i];
                break;
            }
        }
        std::swap(ipnev, ipnod);
        //std::swap(pnev, pnod);
        //std::swap(pnevlen, pnodlen);
        pnev = advec[ipnev];
        pnod = advec[ipnod];
//        pnevlen = adveclen[ipnev]; // not used
//        pnodlen = adveclen[ipnod]; // not used
    }

    adveclen[*iret] = FULL_VEC_LEN;

    ad_free(&ip);
    ad_free(&ipnod);
    ad_free(&ipnev);
}


#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_cos(const TVEC* iv, const TVEC* iret)
#else
void ad_cos(const TVEC* iv, const TVEC* iret)
#endif
{
    //AD<T,nv,nd> ret(v), pnev(v), pnod(v), p(v);
    TVEC ipnev, ipnod, ip;
    ad_alloc(&ipnev);
    ad_alloc(&ipnod);
    ad_alloc(&ip);
    ad_copy(iv, iret);
    ad_copy(iv, &ipnev);
    ad_copy(iv, &ipnod);
    ad_copy(iv, &ip);

    double *v = advec[*iv];
    double *ret = advec[*iret];
    double *pnev = advec[ipnev];
    double *pnod = advec[ipnod];
    double *p = advec[ip];

    // TODO: didn't check c = 0.0
    double s = std::sin(v[0]);
    double c = std::cos(v[0]);

    pnev[0] = pnod[0] = p[0] = 0;

    ret[0] = c;

    size_t plen=adveclen[ip], pnevlen = adveclen[ipnev], pnodlen=adveclen[ipnod];

    for (size_t i = 1; i < plen; ++i) {
        ret[i] *= -s;
    }

    //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
    for (TNVND k = 2; k < gnd+1; ++k) {
        //std::cout << (unsigned int) p.v << std::endl << p << std::endl;
        //std::cout << (unsigned int) pnod.v << std::endl << pnod << std::endl;
        //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
        //AD<T,nv,nd>::mult(p, pnod, pnev);
        ad_mult(&ip, &ipnod, &ipnev);
        //std::cout << (unsigned int) pnev.v << std::endl << pnev << std::endl;
      for (size_t i = 0; i < adveclen[ipnev]; ++i) {
            pnev[i] /= 1.0*k;
            switch(k % 4) {
            case 0:
                ret[i] += c*pnev[i];
                break;
            case 1:
                ret[i] -= s*pnev[i];
                break;
            case 2:
                ret[i] -= c*pnev[i];
                break;
            case 3:
                ret[i] += s*pnev[i];
                break;
            }
        }
        std::swap(ipnev, ipnod);
        pnev = advec[ipnev];
        pnod = advec[ipnod];
        std::swap(pnevlen, pnodlen);
    }

    adveclen[*iret] = FULL_VEC_LEN;

    ad_free(&ip);
    ad_free(&ipnod);
    ad_free(&ipnev);

}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret)
#else
void ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret)
#endif
{
    TNVND* p = base;
    //size_t k = 0;
    unsigned int iexpo = *expo;

    unsigned int *cef = new unsigned int[gnv];
    unsigned int *bv = new unsigned int[gnv];
    unsigned int d = 0, jexp;

    ad_reset(iret);
    advec[*iret][0] = 0;
    adveclen[*iret] = 1;

    // const part was kept
    //if (adveclen[*iv] > 0) {
    //    advec[*iret][0] = advec[*iv][0];
    //    adveclen[*iret] = 1;
    //}

    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        d = 0;
        for (size_t j = 0; j < gnv-1; ++j) {
            //std::cout << " " << c[i];
            cef[j] = *p - *(p+1);
            ++p;
            d += cef[j];
        }
        cef[gnv-1] = *p;
        d += *p;
        ++p;

     #ifdef DEBUG_ALL
        for(size_t j=0; j < gnv; ++j) {
            std::cout << ' ' << cef[j];
        }
        std::cout << "  order: " << d << std::endl;
     #endif

        if (cef[iexpo] <= 0) {
            advec[*iret][i] = 0;
            continue;
        }

        jexp = cef[iexpo];

        cef[iexpo] -= 1;
        --d;

     #ifdef DEBUG_ALL
        std::cout << " --> ";
        for (size_t j = 0; j < gnv; ++j) std::cout << ' ' << cef[j];
        std::cout << "  order: " << d << std::endl;
     #endif

        size_t k = 0;

        //std::cout << "order: " << (unsigned int)d << std::endl;
        for (size_t j = 0; j < gnv; ++j){
            bv[j] = d;
            d -= cef[j];
            k += H[gnv-j][bv[j]];
        }
        //std::cout << std::endl;
        //std::cout << k << std::endl;

        advec[*iret][k] = advec[*iv][i] * 1.0 * jexp;
        if (k >= adveclen[*iret]) adveclen[*iret] = k+1;
     #ifdef DEBUG_ALL
        std::cout << "Set: " << k << ' ' << advec[*iret][k] << "  len: " << adveclen[*iret] << std::endl;
     #endif
    }

    //if (adveclen[*iret] == 0) {
    //    adveclen[*iret] = 1;
    //    advec[*iret] = 0;
    //}

    delete []bv;
    delete []cef;

}

//! similar to ad_derivative, but doesn't multiply the exponent part of x_{expo}
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret)
#else
void ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret)
#endif
{
    TNVND* p = base;
    //size_t k = 0;
    unsigned int iexpo = *expo;

    unsigned int *cef = new unsigned int[gnv];
    unsigned int *bv = new unsigned int[gnv];
    unsigned int d = 0; // , jexp; not used

    ad_reset(iret);
    // const part was kept
    //if (adveclen[*iv] > 0) {
    //    advec[*iret][0] = advec[*iv][0];
    //    adveclen[*iret] = 1;
    //}

    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        d = 0;
        for (size_t j = 0; j < gnv-1; ++j) {
            //std::cout << " " << c[i];
            cef[j] = *p - *(p+1);
            ++p;
            d += cef[j];
        }
        cef[gnv-1] = *p;
        d += *p;
        ++p;

     #ifdef DEBUG_ALL
        for(size_t j=0; j < gnv; ++j) {
            std::cout << ' ' << cef[j];
        }
        std::cout << "  order: " << d << std::endl;
     #endif

        if (cef[iexpo] <= 0) {
            advec[*iret][i] = 0;
            continue;
        }

//        jexp = cef[iexpo]; // not used

        cef[iexpo] -= 1;
        --d;

     #ifdef DEBUG_ALL
        std::cout << " --> ";
        for (size_t j = 0; j < gnv; ++j) std::cout << ' ' << cef[j];
        std::cout << "  order: " << d << std::endl;
     #endif

        size_t k = 0;

        //std::cout << "order: " << (unsigned int)d << std::endl;
        for (size_t j = 0; j < gnv; ++j){
            bv[j] = d;
            d -= cef[j];
            k += H[gnv-j][bv[j]];
        }
        //std::cout << std::endl;
        //std::cout << k << std::endl;

        advec[*iret][k] = advec[*iv][i];
        if (k >= adveclen[*iret]) adveclen[*iret] = k+1;
     #ifdef DEBUG_ALL
        std::cout << "Set: " << k << ' ' << advec[*iret][k] << "  len: " << adveclen[*iret] << std::endl;
     #endif
    }

    delete []bv;
    delete []cef;

}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps)
#else
void ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps)
#endif
{
    double epsv = *eps;
    TNVND* p = base;

    unsigned int *cef = new unsigned int[gnv];
    unsigned int *bv = new unsigned int[gnv];

    ad_reset(iret);

    if (std::abs(*eps) < std::numeric_limits<double>::min()) {
        epsv = std::numeric_limits<double>::min();
    } else {
        epsv = std::abs(*eps);
    }

    unsigned int d = 0, k=0;

    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        if (std::abs(advec[*iv][i]) < epsv) {
            p += gnv;
            continue;
        }

        d = k = 0;
        for (size_t j = 0; j < gnv-1; ++j) {
            //std::cout << " " << c[i];
            cef[j] = *p - *(p+1);
            ++p;
            d += cef[j];
        }
        cef[gnv-1] = *p;
        d += *p;
        ++p;

        for (size_t j = 0; j < *ishift; ++j) {
            k += cef[j];
        }

        if (k > 0) {
            std::cerr << "TPSA: shift a vector, but the first " << *ishift << " components are non-zero" << std::endl;
            std::exit(-1);
        }

        for (size_t j = 0; j < gnv; ++j) {
            if (j+*ishift >= gnv)
                cef[j] = 0;
            else cef[j] = cef[j+*ishift];
        }

     #ifdef DEBUG_ALL
        for(size_t j=0; j < gnv; ++j) {
            std::cout << ' ' << cef[j];
        }
        std::cout << "  order: " << d << std::endl;
     #endif

     #ifdef DEBUG_ALL
        std::cout << " --> ";
        for (size_t j = 0; j < gnv; ++j) std::cout << ' ' << cef[j];
        std::cout << "  order: " << d << std::endl;
     #endif

        k = 0;

        //std::cout << "order: " << (unsigned int)d << std::endl;
        for (size_t j = 0; j < gnv; ++j){
            bv[j] = d;
            d -= cef[j];
            k += H[gnv-j][bv[j]];
        }
        //std::cout << std::endl;
        //std::cout << k << std::endl;

        advec[*iret][k] = advec[*iv][i];
        if (k >= adveclen[*iret]) adveclen[*iret] = k+1;
     #ifdef DEBUG_ALL
        std::cout << "Set: " << k << ' ' << advec[*iret][k] << "  len: " << adveclen[*iret] << std::endl;
     #endif
    }

    delete []bv;
    delete []cef;

}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret)
#else
void ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret)
#endif
{
  (void)nbv; // ldeniau 2011.11.04: avoid g++ warning
  
 #ifdef DEBUG_ALL
    std::cout << "AD substitute " << std::endl;
    std::cout << "  ia[" << adveclen[*iv] << "]:" << *iv
              << " ret:" << *iret << " " << *nbv << std::endl;
    for(size_t i = 0; i < *nbv; ++i) {
        std::cout << " ib  " << i << "  " << ibv[i]
                  << "  len=" << adveclen[ibv[i]] << std::endl;
    }
 #endif

 #ifdef DEBUG_ALL
    std::cout << "v before subst " << std::endl;
    ad_print(iv);
    //std::cout << "ib be new var" << std::endl;
    //for (size_t i = 0; i < *nbv; ++i)
    //    ad_print(&ibv[i]);
 #endif

    ad_reset(iret);
    TVEC it1, it2;
    ad_alloc(&it1);
    ad_alloc(&it2);
    TNVND* pb = base;
    double* pv1 = advec[it1];
    //double* pv2 = advec[it2];
    // loop every entry(monimial)
    for(size_t i = 0; i < adveclen[*iv]; ++i) {
        // only problem of this "== 0" would be(if any) is the speed.
        //if (advec[*iv][i] == 0) {
        if (std::abs(advec[*iv][i]) < std::numeric_limits<double>::min()) {
            pb += gnv;
            continue;
        }

     #ifdef DEBUG_ALL
        std::cout << " base " << i << " : " << *pb << " : ";
     #endif


        // loop over every konob(base vector), nvb == gnv
        ad_reset(&it1);
        ad_reset(&it2);

        pv1[0] = advec[*iv][i];
        adveclen[it1] = 1;
        for(size_t j = 0; j < gnv; ++j) {
            TNVND ord = *(pb+j);
            if (j < gnv-1) ord -= *(pb+j+1);
         #ifdef DEBUG_ALL
            std::cout << ' ' << ord;
         #endif
            if (ord == 0) continue;
            //ad_copy(ibv+j, &it1);
            for (TNVND k = 0; k < ord; ++k) {
                ad_mult(ibv+j, &it1, &it2);
                ad_copy(&it2, &it1);
            }
        }
        //
        ad_add(iret, &it1);
     #ifdef DEBUG_ALL
        std::cout << std::endl;
     #endif
        pb += gnv;
    }
    ad_free(&it1);
    ad_free(&it2);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret)
#else
void ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret)
#endif
{
  (void)iva; (void)nva; (void)iret; (void)nret; // ldeniau 2011.11.04: avoid g++ warnings

 #ifdef DEBUG_ALL
    std::cout << "AD inverse " << std::endl;
    for(size_t i = 0; i < *nva; ++i) {
        std::cout << " iv  " << i << "  " << iva[i] << "  " << iret[i] << std::endl;
        ad_print(&iva[i]);
    }
 #endif

 #ifdef T
    for (TNVND i = 0; i < *nret; ++i) {
        ad_const(&iret[i], 0);
    }

    TVEC* ms = new TVEC[*nret];
    TVEC* ml = new TVEC[*nret];
    TNVND* J = new TNVND[gnv];

    for(size_t i = 0; i < *nret; ++i) {
        ad_alloc(&ms[i]);
        ad_alloc(&ml[i]);
    }

    double vjj;

    for(size_t i = 0; i < *nret; ++i) {
        for (size_t j = 0; j < *nret; ++j) {
            for (size_t k = 0; k < gnv; ++k)
                J[k] = 0;
            J[j] = 1;
            ad_pek(&iva[i], J, &vjj);
            if (std::abs(vjj) > 1e-10) ad_pok(&iva[i], J, gnv, 0);
            //save vjj for inverse.
            if (j!=i){
                ad_pok(&ms[i], J, &gnv, 0.);
                ad_pok(&ml[i], J, &gnv, 0.);
            }else {
                ad_pok(&ms[i], J, &gnv, 1.);
                ad_pok(&ml[i], J, &gnv, 1.);
            }
        }
        //ad_print(&ms[i]);
        //ad_print(&ml[i]);
        ad_mult_const(&iva[i], -1.0);
    }
 #endif

 #ifdef OLD
    // loop gnd order.
    for (size_t iord = 1; iord < gnd; ++iord) {
        for(size_t i = 0; i < *nret; ++i) {
            ad_copy(&vt[i], &vdt[i]);
            // clear the constant part and linear.
            for (size_t k = 0; k < gnv+1; ++k) {
                advec[vdt[i]][k] = 0;
            }
         #ifdef DEBUG_ALL
            std::cout << "vdt  " << i << std::endl;
            ad_print(&vdt[i]);
         #endif
        }
        for(size_t i = 0; i < *nret; ++i) {
            // clear the constant part and linear.
            ad_copy(&vI[i], &vImdt[i]);
            ad_sub(&vImdt[i], &vdt[i]);
        }

        // now vI and vdt is ready
        for (size_t i = 0; i < *nret; ++i) {
            ad_subst(&vt[i], vImdt, nret, &vt1[i]);

         #ifdef DEBUG_ALL
            std::cout << "vt1 after subst " << i << std::endl;
            ad_print(&vt1[i]);
         #endif

        }
        for (size_t i = 0; i < *nret; ++i) {
            ad_copy(&vt1[i], &vt[i]);
        }

        for (size_t i = 0; i < *nret; ++i) {
            ad_subst(&iret[i], vImdt, nret, &vt1[i]);
        }
        for (size_t i = 0; i < *nret; ++i) {
            ad_copy(&vt1[i], &iret[i]);
        }
     #ifdef DEBUG_ALL
        std::cout << "So far the result " << std::endl;
        for (size_t i = 0; i < *nret; ++i)
            ad_print(&iret[i]);
     #endif
    }

    for(size_t i = 0; i < *nret; ++i) {
        ad_free(&vt[i]);
        ad_free(&vt1[i]);
        ad_free(&vI[i]);
        ad_free(&vdt[i]);
        ad_free(&vImdt[i]);
    }
 #endif

 #ifdef T
    for(size_t i = 0; i < *nret; ++i) {
        ad_free(&ml[i]);
        ad_free(&ms[i]);
    }
 #endif
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_nvar(TVEC* n)
#else
void ad_nvar(TVEC* n)
#endif
{
    *n = gnv;
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_length(const TVEC* iv, unsigned int* n)
#else
void ad_length(const TVEC* iv, unsigned int* n)
#endif
{
    *n = adveclen[*iv];
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_read_block(const TVEC* iv,
                                                 double* v, TNVND* J, const unsigned int* N)
#else
void ad_read_block(const TVEC* iv,
                   double* v, TNVND* J, const unsigned int* N)
#endif
{
    //unsigned int ii = *iv;
    if (*N < adveclen[*iv]) {
        for(size_t i = 0; i < *N; ++i) {
            v[i] = 0;
            for (size_t j = 0; j < gnv; ++j) {
                J[i*gnv+j] = 0;
            }
        }
        return;
    }

    TNVND* p = base;
    //os << "iv= " << ii << std::endl;
    //std::cout << "Len: " << *N << ' ' << adveclen[*iv] << std::endl;
    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        v[i] = advec[*iv][i];
        for (size_t j = 0; j < gnv-1; ++j) {
            J[i*gnv+j] = (unsigned int)(*p-*(p+1));
            //std::cout << ' ' << J[i*gnv+j];
            //std::cout << ' ' << *p;
            ++p;
        }
        //std::cout << ' ' << *p << ' ' << &(J[gnv-1][i]) << std::endl;
        J[i*gnv+gnv-1] = *p;
        ++p;
        //std::cout << ' ' << v[i];
    }

    //std::cout << std::endl;

 #ifdef DEBUG_ALL
    std::cout << std::endl;
    for (size_t i = 0; i < adveclen[*iv]*gnv; ++i) {
        std::cout << ' ' << J[i];
        if (i % gnv == gnv-1) std::cout << std::endl;
    }
    //std::cout << "min: " << std::numeric_limits<double>::min() << std::endl;
    //std::cout << "eps: " << std::numeric_limits<double>::epsilon() << std::endl;
 #endif
}

//! blind set, assuming order of J is right.
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_save_block(const TVEC* iv,
                                                 const double* v, const TNVND* J, const unsigned int* N)
#else
void ad_save_block(const TVEC* iv, const double* v, const TNVND* J, const unsigned int* N)
#endif
{
    (void)J; // ldeniau 2011.11.04: avoid g++ warning

    //unsigned int ii = *iv;
    //TNVND* p = base;
    //os << "iv= " << ii << std::endl;
    adveclen[*iv] = *N;

    for (size_t i = 0; i < *N; ++i) {
        advec[*iv][i] = v[i] ;
        //for (size_t j = 0; j < gnv-1; ++j) {
        //    J[i][j] = (*p-*(p+1));
        //    ++p;
        //}
    }

 #ifdef DEBUG_ALL
    //std::cout << "min: " << std::numeric_limits<double>::min() << std::endl;
    //std::cout << "eps: " << std::numeric_limits<double>::epsilon() << std::endl;
 #endif
}

//! Print out ad vector represented by i
#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall print_index(std::ostream& os)
#else
void print_index(std::ostream& os)
#endif
{
    int width_idx=6, width_elem = 3;
    if (FULL_VEC_LEN < 100) width_elem = 3;
    else width_elem = 4;

    os << std::setw(width_idx+width_elem-2) << "Index";
    for (size_t i = 1; i < gnv; ++i) {
        for (size_t j = order_index[i]; j < order_index[i+1]; ++j)
            os << std::setw(width_elem) << j;
        if (i < gnv-1) os << " .";
    }
    os << std::endl;

    os << std::setw(width_idx) << ' ';
    for (size_t i = 0; i < width_elem*order_index[gnv]; ++i) {
        os << '_';
    }
    os << std::endl;

    TNVND ord = 1;
    for (size_t i = 1; i < order_index[gnd]; ++i) {
        if (order_index[ord+1] <= i) ++ord;
        unsigned int M = order_index[gnd-ord+1];
        os << std::setw(width_idx) << i << ' ';

        size_t k = 2;

        for (size_t j = 1; j < M; ++j) {
            if (j == order_index[k]) {
                os << " .";
                ++k;
            }
            os << std::setw(width_elem) << prdidx[i][j];
        }
        std::cout << std::endl;
    }

}

void print_vec(unsigned int ii, std::ostream& os)
{
    //unsigned int ii = *iv;
    TNVND* p = base;
    //os << "iv= " << ii << std::endl;

    std::ios::fmtflags prevflags = os.flags();
    double* v = advec[ii];

    int width_base = 2;
    if (gnd > static_cast<TNVND>(9))  ++width_base;

    os << "          V [" << ii << "]              Base  [ "
       << adveclen[ii] << " / " << FULL_VEC_LEN << " ]" << std::endl
       << "----------------------------------------------" << std::endl;
    for (size_t i = 0; i < adveclen[ii]; ++i) {
        if (std::abs(v[i]) < std::numeric_limits<double>::min()) {
            p += gnv;
            continue;
        }
        os << ' ' << std::setprecision(15)
           << std::scientific << std::setw(15+8) << v[i] << "    ";
        for (size_t j = 0; j < gnv-1; ++j) {
            os << std::setw(width_base) << (unsigned int) (*p-*(p+1));
            ++p;
        }
        os << std::setw(width_base) << (unsigned int)*p++ << std::setw(6) << i << std::endl;
    }
    os << std::endl;

    os.flags(prevflags);

 #ifdef DEBUG_ALL
    //std::cout << "min: " << std::numeric_limits<double>::min() << std::endl;
    //std::cout << "eps: " << std::numeric_limits<double>::epsilon() << std::endl;
 #endif
}


#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_print(const TVEC* iv)
#else
void ad_print(const unsigned int* iv)
#endif
{
    print_vec(*iv, std::cout);
}

#ifdef MSVC_DLL
_declspec(dllexport) void _stdcall ad_print_array(const TVEC* iv, const TVEC* nv)
#else
void ad_print_array(const TVEC* iv, const TVEC* nv)
#endif
{
    std::ostream& os = std::cout;
    //unsigned int ii = *iv;
    //TNVND* p = base;
    //os << "iv= " << ii << std::endl;
 #ifdef V_BY_V
    // vector by vector, each vector one line
    std::ios::fmtflags prevflags = os.flags();
    double* v;
    double x;

    int width_base = 2;
    if (gnd > static_cast<TNVND>(9))  ++width_base;
    os << "    const.  |    linear " << std::endl;
    for (size_t i = 0; i < *nv; ++i) {
        v = advec[iv[i]];
        for (size_t j = 0; j < gnv+1; ++j) {
            if (j < adveclen[iv[i]]) x = v[j];
            else x = 0;

            os << ' ' << std::setprecision(3)
               << std::scientific << std::setw(3+7) << x;
            if (j==0) os << " |";
        }
        os << std::endl;
    }

    os.flags(prevflags);
 #endif

    //unsigned int ii = *iv;
    TNVND* p = base;
    //os << "iv= " << ii << std::endl;
    const char* s = "          ";
    std::ios::fmtflags prevflags = os.flags();
    double* v;
    double xm = 0.0;

    int width_base = 2;
    if (gnd > static_cast<TNVND>(9))  ++width_base;

    std::cout << "# min: " << std::numeric_limits<double>::min() << std::endl;
    std::cout << "# eps: " << std::numeric_limits<double>::epsilon() << std::endl;
    os << "          V              Base  [ "
       << " / " << FULL_VEC_LEN << " ]" << std::endl
       << "----------------------------------------------" << std::endl;
    for (size_t i = 0; i < adveclen[*iv]; ++i) {
        for (size_t j = 0; j < *nv; ++j) {
            v = advec[iv[j]];
            double x = v[i];
            if (std::abs(x) > xm && std::abs(x) < 1)  xm =  std::abs(x);
            if (std::abs(x) < std::numeric_limits<double>::epsilon()){
                x = 0;
                os << ' ' << s;
            } else {
                os << ' ' << std::setprecision(3)
                   << std::scientific << std::setw(3+7) << x;
            }
        }
        os << "   ";
        for (size_t j = 0; j < gnv-1; ++j) {
            os << std::setw(width_base) << (unsigned int) (*p-*(p+1));
            ++p;
        }
        os << std::setw(width_base) << (unsigned int)*p++ << std::setw(6) << i << std::endl;
    }
    os << std::endl;
    os << "# abs(max in [0,1))= " << xm << std::endl;
    os.flags(prevflags);

 #ifdef DEBUG_ALL
 #endif
}

#ifdef MSVC_DLL
#define DLL_PROCESS_DETACH 0
#define DLL_PROCESS_ATTACH 1
extern "C" unsigned long __stdcall DllEntryPoint(void *hDll, unsigned long Reason, void *Reserved)
{
    if (Reason == DLL_PROCESS_ATTACH)
    {
        /*perform DLL initialization tasks here*/
    }
    return (1);
}
#endif
