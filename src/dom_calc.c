/*
 *  Dominance calculation for food web genetic algorithm
 *  Author: Jian Yen
 */

#include<stdio.h>
#include<stdlib.h>
#include<R.h>

void dom_calc(double *fitness, int *n_fit, int *n_pop, int *dom)
{
  // Initialise
  int i, j, k, l;
  int check_k = 1, check_l = 1;

  int num_fit = *n_fit, num_pop = *n_pop;

  int dom_val[num_pop];
  double fitness_vals[num_fit][num_pop];

  // Store dom and fitness in appropriate arrays
  for (i = 0; i < num_pop; i++) {
    dom_val[i] = dom[i];
    for (j = 0; j < num_fit; j++) {
      fitness_vals[j][i] = fitness[j * num_pop + i];
    }
  }

  // Run main dominance calculation
  for (k = 0; k < (num_pop - 1); k++) {
    for (l = k + 1; l < num_pop; l++) {
      check_k = 1;
      check_l = 1;
      for (i = 0; i < num_fit; i++) {
        if (fitness_vals[i][k] < fitness_vals[i][l]) {
          check_l = 0;
        }
        if (fitness_vals[i][l] < fitness_vals[i][k]) {
          check_k = 0;
        }
      }
      if (check_l && !check_k) {
        dom_val[k] += 1;
      }
      if (check_k && !check_l) {
        dom_val[l] += 1;
      }
    }
  }

  //Update and return values
  for (i = 0; i < num_pop; i++) {
    dom[i] = dom_val[i];
    for (j = 0; j < num_fit; j++) {
      fitness[j * num_pop + i] = fitness_vals[j][i];
    }
  }

}