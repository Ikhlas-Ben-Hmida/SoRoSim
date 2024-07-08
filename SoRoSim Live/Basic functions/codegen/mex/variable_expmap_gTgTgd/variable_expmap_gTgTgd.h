/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * variable_expmap_gTgTgd.h
 *
 * Code generation for function 'variable_expmap_gTgTgd'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "variable_expmap_gTgTgd_types.h"

/* Function Declarations */
void variable_expmap_gTgTgd(const real_T Gamma[6], const real_T Gammad[6],
  real_T g[16], real_T Tg[36], real_T Tgd[36]);

/* End of code generation (variable_expmap_gTgTgd.h) */
