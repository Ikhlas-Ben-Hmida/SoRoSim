/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * variable_expmap_gTgTgd_terminate.c
 *
 * Code generation for function 'variable_expmap_gTgTgd_terminate'
 *
 */

/* Include files */
#include "variable_expmap_gTgTgd_terminate.h"
#include "_coder_variable_expmap_gTgTgd_mex.h"
#include "rt_nonfinite.h"
#include "variable_expmap_gTgTgd.h"
#include "variable_expmap_gTgTgd_data.h"

/* Function Definitions */
void variable_expmap_gTgTgd_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void variable_expmap_gTgTgd_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (variable_expmap_gTgTgd_terminate.c) */
