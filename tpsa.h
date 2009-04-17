/*
 * Copyright(C) 2008 by Lingyun Yang
 * lingyun(.dot.]yang@gmail.com
 * http://www.lingyunyang.com
 *
 * Please get permission from Lingyun Yang before you redistribute this file.
 *
 */
//! \brief Automatic Differentiation
//! \file tpsa.h
//! \version $Id: tpsa.h,v 1.4 2009-04-17 17:32:23 frs Exp $
//! \author Lingyun Yang, http://www.lingyunyang.com/

#ifndef WIN32

/* should work unchanged on _win32 using Lahey */
#define ad_print      ad_print_
#define ad_elem       ad_elem_
#define ad_fill_ran   ad_fill_ran_
#define ad_tra        ad_tra_
#define ad_shift      ad_shift_
#define ad_save_block ad_save_block_
#define ad_read_block ad_read_block_
#define ad_nvar       ad_nvar_
#define ad_length     ad_length_
#define ad_derivative ad_derivative_
#define ad_subst      ad_subst_
#define ad_cos        ad_cos_
#define ad_sin        ad_sin_
#define ad_log        ad_log_
#define ad_exp        ad_exp_
#define ad_sqrt       ad_sqrt_
#define ad_abs        ad_abs_
#define ad_div_c      ad_div_c_
#define ad_c_div      ad_c_div_
#define ad_mult_const ad_mult_const_
#define ad_add_const  ad_add_const_
#define ad_div        ad_div_
#define ad_mult       ad_mult_
#define ad_sub        ad_sub_
#define ad_reset      ad_reset_
#define ad_resetvars  ad_resetvars_
#define ad_pok        ad_pok_
#define ad_pek        ad_pek_
#define ad_truncate   ad_truncate_
#define ad_var        ad_var_
#define ad_count      ad_count_
#define ad_const      ad_const_
#define ad_free       ad_free_
#define ad_add        ad_add_
#define ad_copy       ad_copy_
#define ad_clean      ad_clean_
#define ad_alloc      ad_alloc_
#define ad_reserve    ad_reserve_
#define ad_init       ad_init_
#else
#define ad_print      AD_PRINT
#define ad_elem       AD_ELEM
#define ad_fill_ran   AD_FILL_RAN
#define ad_tra        AD_TRA
#define ad_shift      AD_SHIFT
#define ad_save_block AD_SAVE_BLOCK
#define ad_read_block AD_READ_BLOCK
#define ad_nvar       AD_NVAR
#define ad_length     AD_LENGTH
#define ad_derivative AD_DERIVATIVE
#define ad_subst      AD_SUBST
#define ad_cos        AD_COS
#define ad_sin        AD_SIN
#define ad_log        AD_LOG
#define ad_exp        AD_EXP
#define ad_sqrt       AD_SQRT
#define ad_abs        AD_ABS
#define ad_div_c      AD_DIV_C
#define ad_c_div      AD_C_DIV
#define ad_mult_const AD_MULT_CONST
#define ad_add_const  AD_ADD_CONST
#define ad_div        AD_DIV
#define ad_mult       AD_MULT
#define ad_sub        AD_SUB
#define ad_reset      AD_RESET
#define ad_resetvars  AD_RESETVARS_
#define ad_pok        AD_POK
#define ad_pek        AD_PEK
#define ad_truncate   AD_TRUNCATE
#define ad_var        AD_VAR
#define ad_count      AD_COUNT
#define ad_const      AD_CONST
#define ad_free       AD_FREE
#define ad_add        AD_ADD
#define ad_copy       AD_COPY
#define ad_clean      AD_CLEAN
#define ad_alloc      AD_ALLOC
#define ad_reserve    AD_RESERVE
#define ad_init       AD_INIT
#endif


#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#ifndef AD_HH
#define AD_HH
//! Type of order and number of variables.
//typedef unsigned char TNVND;
typedef unsigned int TNVND;
typedef unsigned int TVEC;
#ifdef __cplusplus
extern "C" {
#endif

 #ifdef MSVC_DLL
    _declspec(dllexport) void _stdcall ad_reserve(const TVEC* n);
    _declspec(dllexport) void _stdcall ad_init(const TNVND* nv, const TNVND* nd);
    _declspec(dllexport) void _stdcall ad_resetvars(const TNVND* nv);
    _declspec(dllexport) void _stdcall ad_alloc(TVEC* i);
    _declspec(dllexport) void _stdcall ad_free(const TVEC* i);
    _declspec(dllexport) void _stdcall ad_poolsize(size_t* n);

    _declspec(dllexport) void _stdcall ad_count(TVEC* n);
    _declspec(dllexport) void _stdcall ad_nvar(TVEC* n);
    _declspec(dllexport) void _stdcall ad_length(const TVEC* iv, unsigned int* n);
    _declspec(dllexport) void _stdcall ad_copy(const TVEC* i, const TVEC* j);
    _declspec(dllexport) void _stdcall ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x);
    _declspec(dllexport) void _stdcall ad_pek(const TVEC* ivec, int* c, size_t* n, double* x);
    _declspec(dllexport) void _stdcall ad_pok(const TVEC* ivec, int* c, size_t* n, double* x);
    _declspec(dllexport) void _stdcall ad_var(const TVEC* ii, const double* x, unsigned int* iv);
    _declspec(dllexport) void _stdcall ad_abs(const TVEC* iv, double* r);
    _declspec(dllexport) void _stdcall ad_truncate(const TVEC* iv, const TNVND* d);

    _declspec(dllexport) void _stdcall ad_clean(const TVEC* iv, const double* eps);
    _declspec(dllexport) void _stdcall ad_reset(const TVEC* iv);
    _declspec(dllexport) void _stdcall ad_const(const TVEC* ii, const double* r);
    _declspec(dllexport) void _stdcall ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm);

    _declspec(dllexport) void _stdcall ad_add(const TVEC* i, const TVEC* j);
    _declspec(dllexport) void _stdcall ad_add_const(const TVEC* i, double *r);

    _declspec(dllexport) void _stdcall ad_sub(const TVEC* i, const TVEC* j);
    //! internal multiplication, dst should be different from lhs and rhs.
    _declspec(dllexport) void _stdcall ad_mult(const TVEC* ivlhs, const TVEC* ivrhs, TVEC* ivdst);
    _declspec(dllexport) void _stdcall ad_mult_const(const TVEC* iv, double* c);

    _declspec(dllexport) void _stdcall ad_div_c(const TVEC* iv, const double* c);
    _declspec(dllexport) void _stdcall ad_c_div(const TVEC* iv, const double* c, TVEC* ivret);
    _declspec(dllexport) void _stdcall ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst);

    _declspec(dllexport) void _stdcall ad_sqrt(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_exp(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_log(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_sin(const TVEC* iv, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_cos(const TVEC* iv, const TVEC* iret);

    _declspec(dllexport) void _stdcall ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    _declspec(dllexport) void _stdcall ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps);

    _declspec(dllexport) void _stdcall ad_read_block(const TVEC* iv, double* v, TNVND* J, const unsigned int* N);
    _declspec(dllexport) void _stdcall ad_save_block(const TVEC* iv, const double* v, const TNVND* J, const unsigned int* N);

    _declspec(dllexport) void _stdcall ad_rev(const TVEC* iv);

    _declspec(dllexport) void _stdcall ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret);
    //void ad_inverse(const TVEC* iv, const TNVND* nbv, const TVEC* iret, const TNVND* nret);
    _declspec(dllexport) void _stdcall ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret);

    _declspec(dllexport) void _stdcall print_index(std::ostream& os);
    _declspec(dllexport) void _stdcall ad_print(const TVEC* iv);
    _declspec(dllexport) void _stdcall ad_print_array(const TVEC* iv, const TVEC* nv);
 #else
    void ad_reserve(const TVEC* n);
    void ad_init(const TNVND* nv, const TNVND* nd);
    void ad_resetvars(const TNVND* nv);
    void ad_alloc(TVEC* i);
    void ad_free(const TVEC* i);
    void ad_poolsize(size_t* n);

    void ad_count(TVEC* n);
    void ad_nvar(TVEC* n);
    void ad_length(const TVEC* iv, unsigned int* n);
    void ad_copy(const TVEC* i, const TVEC* j);
    void ad_elem(const TVEC* ivec, unsigned int* idx, unsigned int* c, double* x);
    void ad_pek(const TVEC* ivec, int* c, size_t* n, double* x);
    void ad_pok(const TVEC* ivec, int* c, size_t* n, double* x);
    void ad_var(const TVEC* ii, const double* x, unsigned int* iv);
    void ad_abs(const TVEC* iv, double* r);
    void ad_truncate(const TVEC* iv, const TNVND* d);

    void ad_clean(const TVEC* iv, const double* eps);
    void ad_reset(const TVEC* iv);
    void ad_const(const TVEC* ii, const double* r);
    void ad_fill_ran(const TVEC* iv, const double* ratio, const double* xm);

    void ad_add(const TVEC* i, const TVEC* j);
    void ad_add_const(const TVEC* i, double *r);

    void ad_sub(const TVEC* i, const TVEC* j);
    //! internal multiplication, dst should be different from lhs and rhs.
    void ad_mult(const TVEC* ivlhs, const TVEC* ivrhs, TVEC* ivdst);
    void ad_mult_const(const TVEC* iv, double* c);

    void ad_div_c(const TVEC* iv, const double* c);
    void ad_c_div(const TVEC* iv, const double* c, TVEC* ivret);
    void ad_div(const TVEC* ilhs, const TVEC* irhs, TVEC* idst);

    void ad_sqrt(const TVEC* iv, const TVEC* iret);
    void ad_exp(const TVEC* iv, const TVEC* iret);
    void ad_log(const TVEC* iv, const TVEC* iret);
    void ad_sin(const TVEC* iv, const TVEC* iret);
    void ad_cos(const TVEC* iv, const TVEC* iret);

    void ad_derivative(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    void ad_tra(const TVEC* iv, unsigned int* expo, const TVEC* iret);
    void ad_shift(const TVEC* iv, unsigned int* ishift, const TVEC* iret, const double* eps);

    void ad_read_block(const TVEC* iv, double* v, TNVND* J, const unsigned int* N);
    void ad_save_block(const TVEC* iv, const double* v, const TNVND* J, const unsigned int* N);

    void ad_rev(const TVEC* iv);

    void ad_subst(const TVEC* iv, const TVEC* ibv, const TNVND* nbv, const TVEC* iret);
    void ad_inverse(const TVEC* iva, const TNVND* nva, const TVEC* iret, const TNVND* nret);

    void print_index(std::ostream& os);
    void ad_print(const TVEC* iv);
    void ad_print_array(const TVEC* iv, const TVEC* nv);
 #endif

#ifdef __cplusplus
}
#endif

#endif
