/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rt_nonfinite.c
 *
 * Code generation for function 'variable_expmap_gTgTgd'
 *
 */

#include "rt_nonfinite.h"
#include "matrix.h"

real_T mex_rtInf;
real_T mex_rtMinusInf;
real_T mex_rtNaN;
real32_T mex_rtInfF;
real32_T mex_rtMinusInfF;
real32_T mex_rtNaNF;
void mex_InitInfAndNan(void)
{
  mex_rtInf = mxGetInf();
  mex_rtMinusInf = -rtInf;
  mex_rtInfF = (real32_T)rtInf;
  mex_rtMinusInfF = -rtInfF;
  mex_rtNaN = mxGetNaN();
  mex_rtNaNF = (real32_T)rtNaN;
}

/* End of code generation (rt_nonfinite.c) */
