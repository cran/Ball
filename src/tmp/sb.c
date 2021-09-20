#include "bf_c.h"
#include "stdlib.h"
#include "stdio.h"
#include "utilities.h"
#include "R.h"

double sbd(double **xx, double **xy, double **yx, double **yy, int n1, int n1_total, int n2, int n2_total) {
  void *bf_x = BF_C_new(xx, n1, n1_total);
  void *bf_y = BF_C_new(yy, n2, n2_total);
  // Rprintf("initial pass\n");
  
  double **pxx = alloc_matrix(n1, n1_total);
  double **pxy = alloc_matrix(n1, n1_total);
  double **pyx = alloc_matrix(n2, n2_total);
  double **pyy = alloc_matrix(n2, n2_total);

  BF_C_get_fitted(bf_x, pxx);
  BF_C_get_fitted(bf_y, pyy);
  // Rprintf("get fitted pass\n");

  BF_C_predict(bf_x, pxy, xy, n2_total);
  // Rprintf("\n predict bf_x pass\n");
  BF_C_predict(bf_y, pyx, yx, n1_total);
  // Rprintf("\n predict bf_y pass\n");
  // Rprintf("get predict pass\n");
  
  double tmp;
  double denominator; 
  double A = 0.0, C = 0.0, bd;

  // Rprintf("Pxx matrix: \n");
  // for (int i = 0; i < n1; i++) {
  //   for (int j = 0; j < n1_total; j++) {
  //     Rprintf("%f ", pxx[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  // Rprintf("Pxy matrix: \n");
  // for (int i = 0; i < n1; i++) {
  //   for (int j = 0; j < n1_total; j++) {
  //     Rprintf("%f ", pxy[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  // Rprintf("Pyy matrix: \n");
  // for (int i = 0; i < n2; i++) {
  //   for (int j = 0; j < n2_total; j++) {
  //     Rprintf("%f ", pyy[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  // Rprintf("Pyx matrix: \n");
  // for (int i = 0; i < n2; i++) {
  //   for (int j = 0; j < n1_total; j++) {
  //     Rprintf("%f ", pyx[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n1_total; j++) {
      tmp = pxx[i][j] - pxy[i][j];
      A += tmp * tmp;
    }
  }
  denominator = 1.0 / ((double) n1 * n1_total);
  A *= denominator;
  
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < n2_total; j++) {
      tmp = pyy[i][j] - pyx[i][j];
      C += tmp * tmp;
    }
  }
  denominator = 1.0 / ((double) n2 * n2_total);
  C *= denominator;

  bd = A + C;
  // Rprintf("A: %f; C: %f; BD: %f\n", A, C, bd);

  free_matrix(pxx, n1, n1_total);
  free_matrix(pxy, n1, n1_total);
  free_matrix(pyx, n2, n2_total);
  free_matrix(pyy, n2, n2_total);
  // Rprintf("probability matrix \n");
  
  BF_C_free_BF(bf_x);
  BF_C_free_BF(bf_y);
  // Rprintf("free BF's pointer \n");

  BF_C_delete(bf_x);
  BF_C_delete(bf_y);
  // Rprintf("destruct bf_x and bf_y\n");

  return bd;
}

void sbd_C(double *sbd_value, double *x, int *n1, int *n1_total, int *n2, int *n2_total) {
  int num = *n1_total + *n2_total;
  int fix_center_num = *n1 + *n2;

  double **xx = alloc_matrix(*n1, *n1_total);
  double **xy = alloc_matrix(*n1, *n2_total);
  double **yx = alloc_matrix(*n2, *n1_total);
  double **yy = alloc_matrix(*n2, *n2_total);

  for (int i = 0; i < *n1; i++) {
    for (int j = 0; j < *n1_total; j++) {
      xx[i][j] = x[i*num + j];
    }
  }

  for (int i = 0; i < *n1; i++) {
    for (int j = *n1_total; j < num; j++) {
      xy[i][j - *n1_total] = x[i*num + j];
    }
  }

  for (int i = *n1; i < fix_center_num; i++) {
    for (int j = 0; j < *n1_total; j++) {
      yx[i - *n1][j] = x[i*num + j];
    }
  }

  for (int i = *n1; i < fix_center_num; i++) {
    for (int j = *n1_total; j < num; j++) {
      yy[i - *n1][j - *n1_total] = x[i*num + j];
    }
  }

  // Rprintf("xx matrix: \n");
  // for (int i = 0; i < *n1; i++) {
  //   for (int j = 0; j < *n1_total; j++) {
  //     Rprintf("%f ", xx[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  // Rprintf("xy matrix: \n");
  // for (int i = 0; i < *n1; i++) {
  //   for (int j = 0; j < *n2_total; j++) {
  //     Rprintf("%f ", xy[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  *sbd_value = sbd(xx, xy, yx, yy, *n1, *n1_total, *n2, *n2_total);
  // Rprintf("compute sbd value successful.\n");

  free_matrix(xx, *n1, *n1_total);
  free_matrix(xy, *n1, *n2_total);
  free_matrix(yx, *n2, *n1_total);
  free_matrix(yy, *n2, *n2_total);

  return;
}

// double sbd(double **xx, double **yy, double **xy, double **yx) {
//   BF bf_x = BF_C_new(xx);
//   BF bf_y = BF(yy);
  
//   double **pxx = bf_x.predict();
//   double **pxy = bf_x.predict(xy);
  
//   double **pyy = bf_y.predict();
//   double **pyx = bf_y.predict(yx);
  
//   unsigned int n1, n2;
//   double bd = 0.0;
  
//   n1 = sizeof(pxx);
//   n2 = sizeof(pxx[0]);
//   double tmp;
//   for (unsigned int i = 0; i < n1; i++) {
//     for (unsigned int j = 0; j < n2; j++) {
//       tmp = pxx[i][j] - pxy[i][j];
//       bd += tmp * tmp;
//     }
//   }
  
//   n1 = sizeof(pyy);
//   n2 = sizeof(pyy[0]);
//   for (unsigned int i = 0; i < n1; i++) {
//     for (unsigned int j = 0; j < n2; j++) {
//       tmp = pyy[i][j] - pyx[i][j];
//       bd += tmp * tmp;
//     }
//   }
  
//   return bd;
// }
