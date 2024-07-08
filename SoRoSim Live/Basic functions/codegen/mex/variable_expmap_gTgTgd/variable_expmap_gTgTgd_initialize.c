/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * variable_expmap_gTgTgd_initialize.c
 *
 * Code generation for function 'variable_expmap_gTgTgd_initialize'
 *
 */

/* Include files */
#include "variable_expmap_gTgTgd_initialize.h"
#include "_coder_variable_expmap_gTgTgd_mex.h"
#include "rt_nonfinite.h"
#include "variable_expmap_gTgTgd.h"
#include "variable_expmap_gTgTgd_data.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;

/* Function Definitions */
void variable_expmap_gTgTgd_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (variable_expmap_gTgTgd_initialize.c) */
