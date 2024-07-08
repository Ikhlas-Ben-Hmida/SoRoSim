/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rt_nonfinite.h
 *
 * Code generation for function 'variable_expmap_gTgTgd'
 *
 */

#pragma once

#include "tmwtypes.h"

extern real_T mex_rtInf;
extern real_T mex_rtMinusInf;
extern real_T mex_rtNaN;
extern real32_T mex_rtInfF;
extern real32_T mex_rtMinusInfF;
extern real32_T mex_rtNaNF;

#define rtInf                          mex_rtInf
#define rtMinusInf                     mex_rtMinusInf
#define rtNaN                          mex_rtNaN
#define rtInfF                         mex_rtInfF
#define rtMinusInfF                    mex_rtMinusInfF
#define rtNaNF                         mex_rtNaNF
#define rtIsNaN(X)                     mxIsNaN(X)
#define rtIsInf(X)                     mxIsInf(X)
#define rtIsNaNF(X)                    mxIsNaN(X)
#define rtIsInfF(X)                    mxIsInf(X)

extern void mex_InitInfAndNan(void);

/* End of code generation (rt_nonfinite.h) */
