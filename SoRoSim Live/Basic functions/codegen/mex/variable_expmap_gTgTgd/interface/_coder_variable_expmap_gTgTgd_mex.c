/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_variable_expmap_gTgTgd_mex.c
 *
 * Code generation for function '_coder_variable_expmap_gTgTgd_mex'
 *
 */

/* Include files */
#include "_coder_variable_expmap_gTgTgd_mex.h"
#include "_coder_variable_expmap_gTgTgd_api.h"
#include "variable_expmap_gTgTgd.h"
#include "variable_expmap_gTgTgd_data.h"
#include "variable_expmap_gTgTgd_initialize.h"
#include "variable_expmap_gTgTgd_terminate.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void c_variable_expmap_gTgTgd_mexFun(int32_T nlhs, mxArray
  *plhs[3], int32_T nrhs, const mxArray *prhs[2]);

/* Function Definitions */
void c_variable_expmap_gTgTgd_mexFun(int32_T nlhs, mxArray *plhs[3], int32_T
  nrhs, const mxArray *prhs[2])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        22, "variable_expmap_gTgTgd");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "variable_expmap_gTgTgd");
  }

  /* Call the function. */
  variable_expmap_gTgTgd_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&variable_expmap_gTgTgd_atexit);

  /* Module initialization. */
  variable_expmap_gTgTgd_initialize();

  /* Dispatch the entry-point. */
  c_variable_expmap_gTgTgd_mexFun(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  variable_expmap_gTgTgd_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_variable_expmap_gTgTgd_mex.c) */
