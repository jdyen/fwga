/*
 *  Dominance calculation for food web genetic algorithm
 *  Author: Jian Yen
 */

#include<stdio.h>
#include<stdlib.h>
#include<R.h>

void fitness_se(double *fw_in, int *n_sp, double *sec_e)
{
  // Initialise
  int i, j, k;

  int n_species = *n_sp;

  double sec_ext = *sec_e;
  double fw[n_species][n_species];
  double fw_store[n_species][n_species];
  double col_sum = 0;

  // Store fw and prim_prod in appropriate arrays
  for (i = 0; i < n_species; i++) {
    for (j = 0; j < n_species; j++) {
      fw[i][j] = fw_in[i * n_species + j];
    }
    fw[i][i] = 0;
  }

  // Run secondary extinction calculation
  for (i = 0; i < n_species; i++) {
    for (j = 0; j < n_species; j++) {
      fw_store[i][j] = fw[i][j];
      fw[i][j] = 0;
      col_sum = 0;
      for (k = 0; k < n_species; k++) {
        col_sum += fw[k][j];
      }
      if ((col_sum == 0) && (j != i)) {
        sec_ext += 1;
      }
      fw[i][j] = fw_store[i][j];
    }
  }

  //Update and return values
  for (i = 0; i < n_species; i++) {
    for (j = 0; j < n_species; j++) {
      fw_in[i * n_species + j] = fw[i][j];
    }
  }
  *sec_e = sec_ext;

}