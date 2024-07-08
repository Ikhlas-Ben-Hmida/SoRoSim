/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_variable_expmap_gTgTgd_api.c
 *
 * Code generation for function '_coder_variable_expmap_gTgTgd_api'
 *
 */

/* Include files */
#include "_coder_variable_expmap_gTgTgd_api.h"
#include "rt_nonfinite.h"
#include "variable_expmap_gTgTgd.h"
#include "variable_expmap_gTgTgd_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6];
static const mxArray *b_emlrt_marshallOut(const real_T u[36]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *Gamma,
  const char_T *identifier))[6];
static const mxArray *emlrt_marshallOut(const real_T u[16]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6]
{
  real_T (*y)[6];
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(const real_T u[36])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 6, 6 };

  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6]
{
  real_T (*ret)[6];
  static const int32_T dims[1] = { 6 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[6])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *Gamma,
  const char_T *identifier))[6]
{
  real_T (*y)[6];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(Gamma), &thisId);
  emlrtDestroyArray(&Gamma);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u[16])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 4, 4 };

  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

void variable_expmap_gTgTgd_api(const mxArray * const prhs[2], int32_T nlhs,
  const mxArray *plhs[3])
{
  real_T (*g)[16];
  real_T (*Tg)[36];
  real_T (*Tgd)[36];
  real_T (*Gamma)[6];
  real_T (*Gammad)[6];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  g = (real_T (*)[16])mxMalloc(sizeof(real_T [16]));
  Tg = (real_T (*)[36])mxMalloc(sizeof(real_T [36]));
  Tgd = (real_T (*)[36])mxMalloc(sizeof(real_T [36]));

  /* Marshall function inputs */
  Gamma = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "Gamma");
  Gammad = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "Gammad");

  /* Invoke the target function */
  variable_expmap_gTgTgd(*Gamma, *Gammad, *g, *Tg, *Tgd);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*g);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*Tg);
  }

  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(*Tgd);
  }
}

/* End of code generation (_coder_variable_expmap_gTgTgd_api.c) */
