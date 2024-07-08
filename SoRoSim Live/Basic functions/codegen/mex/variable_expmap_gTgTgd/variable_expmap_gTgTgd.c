/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * variable_expmap_gTgTgd.c
 *
 * Code generation for function 'variable_expmap_gTgTgd'
 *
 */

/* Include files */
#include "variable_expmap_gTgTgd.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void variable_expmap_gTgTgd(const real_T Gamma[6], const real_T Gammad[6],
  real_T g[16], real_T Tg[36], real_T Tgd[36])
{
  real_T scale;
  real_T absxk;
  real_T t;
  real_T theta;
  real_T thetad;
  real_T adjGammad[36];
  int32_T i;
  int32_T i1;
  static const int8_T iv[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
  };

  int32_T adjGamma_tmp;
  static const int8_T iv1[36] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  real_T Gammahatp2[16];
  real_T adjGammad2[36];
  real_T adjGamma[36];
  real_T d;
  real_T adjGammap2[36];
  real_T tp2;
  real_T adjGammad3[36];
  real_T tp3;
  real_T d1;
  real_T adjGammap3[36];
  real_T tp4;
  real_T tp5;
  real_T sintheta;
  int32_T i2;
  int32_T i3;
  real_T t2;
  real_T t5;
  real_T adjGammap4[36];
  real_T t6;
  real_T t7;
  real_T t8;
  real_T a;
  real_T b_a;
  real_T b_Gammahatp2[16];
  real_T b_adjGammap3[36];
  scale = 3.3121686421112381E-170;
  absxk = muDoubleScalarAbs(Gamma[0]);
  if (absxk > 3.3121686421112381E-170) {
    theta = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    theta = t * t;
  }

  absxk = muDoubleScalarAbs(Gamma[1]);
  if (absxk > scale) {
    t = scale / absxk;
    theta = theta * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    theta += t * t;
  }

  absxk = muDoubleScalarAbs(Gamma[2]);
  if (absxk > scale) {
    t = scale / absxk;
    theta = theta * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    theta += t * t;
  }

  theta = scale * muDoubleScalarSqrt(theta);
  thetad = ((Gammad[0] * Gamma[0] + Gammad[1] * Gamma[1]) + Gammad[2] * Gamma[2])
    / theta;

  /*  optimized on 30.05.2022 */
  g[0] = 0.0;
  g[4] = -Gamma[2];
  g[8] = Gamma[1];
  g[12] = Gamma[3];
  g[1] = Gamma[2];
  g[5] = 0.0;
  g[9] = -Gamma[0];
  g[13] = Gamma[4];
  g[2] = -Gamma[1];
  g[6] = Gamma[0];
  g[10] = 0.0;
  g[14] = Gamma[5];
  g[3] = 0.0;
  g[7] = 0.0;
  g[11] = 0.0;
  g[15] = 0.0;

  /*  optimized on 30.05.2022 */
  Tgd[0] = 0.0;
  Tgd[6] = -Gamma[2];
  Tgd[12] = Gamma[1];
  Tgd[18] = 0.0;
  Tgd[24] = 0.0;
  Tgd[30] = 0.0;
  Tgd[1] = Gamma[2];
  Tgd[7] = 0.0;
  Tgd[13] = -Gamma[0];
  Tgd[19] = 0.0;
  Tgd[25] = 0.0;
  Tgd[31] = 0.0;
  Tgd[2] = -Gamma[1];
  Tgd[8] = Gamma[0];
  Tgd[14] = 0.0;
  Tgd[20] = 0.0;
  Tgd[26] = 0.0;
  Tgd[32] = 0.0;
  Tgd[3] = 0.0;
  Tgd[9] = -Gamma[5];
  Tgd[15] = Gamma[4];
  Tgd[21] = 0.0;
  Tgd[27] = -Gamma[2];
  Tgd[33] = Gamma[1];
  Tgd[4] = Gamma[5];
  Tgd[10] = 0.0;
  Tgd[16] = -Gamma[3];
  Tgd[22] = Gamma[2];
  Tgd[28] = 0.0;
  Tgd[34] = -Gamma[0];
  Tgd[5] = -Gamma[4];
  Tgd[11] = Gamma[3];
  Tgd[17] = 0.0;
  Tgd[23] = -Gamma[1];
  Tgd[29] = Gamma[0];
  Tgd[35] = 0.0;

  /*  optimized on 30.05.2022 */
  adjGammad[0] = 0.0;
  adjGammad[6] = -Gammad[2];
  adjGammad[12] = Gammad[1];
  adjGammad[18] = 0.0;
  adjGammad[24] = 0.0;
  adjGammad[30] = 0.0;
  adjGammad[1] = Gammad[2];
  adjGammad[7] = 0.0;
  adjGammad[13] = -Gammad[0];
  adjGammad[19] = 0.0;
  adjGammad[25] = 0.0;
  adjGammad[31] = 0.0;
  adjGammad[2] = -Gammad[1];
  adjGammad[8] = Gammad[0];
  adjGammad[14] = 0.0;
  adjGammad[20] = 0.0;
  adjGammad[26] = 0.0;
  adjGammad[32] = 0.0;
  adjGammad[3] = 0.0;
  adjGammad[9] = -Gammad[5];
  adjGammad[15] = Gammad[4];
  adjGammad[21] = 0.0;
  adjGammad[27] = -Gammad[2];
  adjGammad[33] = Gammad[1];
  adjGammad[4] = Gammad[5];
  adjGammad[10] = 0.0;
  adjGammad[16] = -Gammad[3];
  adjGammad[22] = Gammad[2];
  adjGammad[28] = 0.0;
  adjGammad[34] = -Gammad[0];
  adjGammad[5] = -Gammad[4];
  adjGammad[11] = Gammad[3];
  adjGammad[17] = 0.0;
  adjGammad[23] = -Gammad[1];
  adjGammad[29] = Gammad[0];
  adjGammad[35] = 0.0;
  if (theta <= 1.0E-6) {
    for (i = 0; i < 16; i++) {
      g[i] += (real_T)iv[i];
    }

    for (i = 0; i < 36; i++) {
      Tg[i] = (real_T)iv1[i] + 0.5 * Tgd[i];
      Tgd[i] = 0.5 * adjGammad[i];
    }
  } else {
    for (i = 0; i < 4; i++) {
      for (i1 = 0; i1 < 4; i1++) {
        adjGamma_tmp = i1 << 2;
        Gammahatp2[i + adjGamma_tmp] = ((g[i] * g[adjGamma_tmp] + g[i + 4] *
          g[adjGamma_tmp + 1]) + g[i + 8] * g[adjGamma_tmp + 2]) + g[i + 12] *
          g[adjGamma_tmp + 3];
      }
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        for (adjGamma_tmp = 0; adjGamma_tmp < 6; adjGamma_tmp++) {
          d += Tgd[i + 6 * adjGamma_tmp] * Tgd[adjGamma_tmp + 6 * i1];
        }

        adjGammap2[i + 6 * i1] = d;
      }

      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        for (adjGamma_tmp = 0; adjGamma_tmp < 6; adjGamma_tmp++) {
          d += adjGammap2[i + 6 * adjGamma_tmp] * Tgd[adjGamma_tmp + 6 * i1];
        }

        adjGammap3[i + 6 * i1] = d;
      }

      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        d1 = 0.0;
        scale = 0.0;
        for (adjGamma_tmp = 0; adjGamma_tmp < 6; adjGamma_tmp++) {
          i2 = adjGamma_tmp + 6 * i1;
          i3 = i + 6 * adjGamma_tmp;
          d += adjGammap3[i3] * Tgd[i2];
          d1 += adjGammad[i3] * Tgd[i2];
          scale += Tgd[i3] * adjGammad[i2];
        }

        adjGamma_tmp = i + 6 * i1;
        adjGamma[adjGamma_tmp] = scale;
        adjGammad2[adjGamma_tmp] = d1;
        adjGammap4[adjGamma_tmp] = d;
      }
    }

    for (i = 0; i < 36; i++) {
      adjGammad2[i] += adjGamma[i];
    }

    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        d1 = 0.0;
        for (adjGamma_tmp = 0; adjGamma_tmp < 6; adjGamma_tmp++) {
          i2 = i + 6 * adjGamma_tmp;
          i3 = adjGamma_tmp + 6 * i1;
          d += adjGammad2[i2] * Tgd[i3];
          d1 += adjGammap2[i2] * adjGammad[i3];
        }

        adjGamma_tmp = i + 6 * i1;
        adjGamma[adjGamma_tmp] = d1;
        adjGammad3[adjGamma_tmp] = d;
      }
    }

    for (i = 0; i < 36; i++) {
      adjGammad3[i] += adjGamma[i];
    }

    tp2 = theta * theta;
    tp3 = tp2 * theta;
    tp4 = tp3 * theta;
    tp5 = tp4 * theta;
    sintheta = muDoubleScalarSin(theta);
    absxk = muDoubleScalarCos(theta);
    t = theta * sintheta;
    t2 = theta * absxk;
    t5 = ((4.0 - 4.0 * absxk) - t) / (2.0 * tp2);
    t6 = ((4.0 * theta - 5.0 * sintheta) + t2) / (2.0 * tp3);
    t7 = ((2.0 - 2.0 * absxk) - t) / (2.0 * tp4);
    t8 = ((2.0 * theta - 3.0 * sintheta) + t2) / (2.0 * tp5);
    a = (1.0 - absxk) / tp2;
    b_a = (theta - sintheta) / tp3;
    for (i = 0; i < 4; i++) {
      d = Gammahatp2[i + 4];
      d1 = Gammahatp2[i + 8];
      scale = Gammahatp2[i + 12];
      for (i1 = 0; i1 < 4; i1++) {
        adjGamma_tmp = i1 << 2;
        b_Gammahatp2[i + adjGamma_tmp] = ((Gammahatp2[i] * g[adjGamma_tmp] + d *
          g[adjGamma_tmp + 1]) + d1 * g[adjGamma_tmp + 2]) + scale *
          g[adjGamma_tmp + 3];
      }
    }

    for (i = 0; i < 16; i++) {
      g[i] = (((real_T)iv[i] + g[i]) + a * Gammahatp2[i]) + b_a * b_Gammahatp2[i];
    }

    for (i = 0; i < 36; i++) {
      Tg[i] = ((((real_T)iv1[i] + t5 * Tgd[i]) + t6 * adjGammap2[i]) + t7 *
               adjGammap3[i]) + t8 * adjGammap4[i];
    }

    absxk = thetad * (((8.0 - tp2) * absxk + -8.0) + 5.0 * t);
    a = absxk / (2.0 * tp3);
    scale = thetad * ((-8.0 * theta + (15.0 - tp2) * sintheta) - 7.0 * t2);
    b_a = scale / (2.0 * tp4);
    absxk /= 2.0 * tp5;
    scale /= 2.0 * (tp5 * theta);
    for (i = 0; i < 6; i++) {
      for (i1 = 0; i1 < 6; i1++) {
        d = 0.0;
        d1 = 0.0;
        for (adjGamma_tmp = 0; adjGamma_tmp < 6; adjGamma_tmp++) {
          i2 = i + 6 * adjGamma_tmp;
          i3 = adjGamma_tmp + 6 * i1;
          d += adjGammad3[i2] * Tgd[i3];
          d1 += adjGammap3[i2] * adjGammad[i3];
        }

        adjGamma_tmp = i + 6 * i1;
        b_adjGammap3[adjGamma_tmp] = d1;
        adjGamma[adjGamma_tmp] = d;
      }
    }

    for (i = 0; i < 36; i++) {
      Tgd[i] = ((((((a * Tgd[i] + t5 * adjGammad[i]) + b_a * adjGammap2[i]) + t6
                   * adjGammad2[i]) + absxk * adjGammap3[i]) + t7 * adjGammad3[i])
                + scale * adjGammap4[i]) + t8 * (adjGamma[i] + b_adjGammap3[i]);
    }
  }
}

/* End of code generation (variable_expmap_gTgTgd.c) */
