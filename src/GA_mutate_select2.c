/*
 *  Dominance calculation for food web genetic algorithm
 *  Author: Jian Yen
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void GA_mutate_select2(double *pops, int *n_species, int *n_ones, int *n_zero,
                       double *damp_val)
{
  // Initialise
  int i, j;
  int num_species = *n_species;
  int ones_count = *n_ones, zero_count = *n_zero;
  double damper = *damp_val;

  int change_count, change_range, num_changed, num_tested, plus_minus, num_to_change;
  double a_select, select_prob, increment;
  
  double new_pop[num_species * num_species], store_pop[num_species * num_species];
  for (i = 0; i < (num_species * num_species); i++) {
    new_pop[i] = pops[i];
    store_pop[i] = pops[i];
  }
    
  // choose change_count links to change in the candidate pop
  if (zero_count < ones_count) {
    change_range = zero_count;
  } else {
    change_range = ones_count;
  }
  double change_probs[change_range];
  change_probs[0] = 1.0 / change_range;
  for (i = 1; i < change_range; i++) {
    change_probs[i] = change_probs[i - 1] + 1.0 / change_range; 
  }
  GetRNGstate();
  a_select = runif(0.0, 1.0);
  PutRNGstate();
  change_count = 0;
  while (a_select > change_probs[change_count]) {
    change_count += 1;
  }

  // select change_count zero links to make positive
  num_changed = 0;
  num_tested = 0;
  i = 0;
  j = 0;
  while (num_changed < change_count) {
    if (new_pop[i * num_species + j] == 0.0) {
      select_prob = (double)(change_count - num_changed) / (double)(zero_count - num_tested);
      GetRNGstate();
      a_select = runif(0.0, 1.0);
      PutRNGstate();
      if (a_select < select_prob) {
        GetRNGstate();
        new_pop[i * num_species + j] = runif(0.01, 1.0);
        PutRNGstate();
        num_changed += 1;
      }
      num_tested += 1;
    }
    if ((i < num_species) && (j < (num_species - 1))) {
      j += 1;
    } else {
      if (j == (num_species - 1)) {
        j = 0;
        i += 1;
      }
    }
  }

  // select change_count positive links to make zero
  num_changed = 0;
  num_tested = 0;
  i = 0;
  j = 0;
  while (num_changed < change_count) {
    if (store_pop[i * num_species + j] != 0.0) {
      select_prob = (double)(change_count - num_changed) / (double)(ones_count - num_tested);
      GetRNGstate();
      a_select = runif(0.0, 1.0);
      PutRNGstate();
      if (a_select < select_prob) {
        new_pop[i * num_species + j] = 0;
        num_changed += 1;
      }
      num_tested += 1;
    }
    if ((i < num_species) && (j < (num_species - 1))) {
      j += 1;
    } else {
      if (j == (num_species - 1)) {
        j = 0;
        i += 1;
      }
    }
  }
    
  // select randomly from positive links and add or subtract the damper value
  num_changed = 0;
  num_tested = 0;
  num_to_change = (int) floor(ones_count / 3.0);
  i = 0;
  j = 0;
  while (num_changed < num_to_change) {
    if (new_pop[i * num_species + j] != 0.0) {
      select_prob = (double)(num_to_change - num_changed) / (double)(ones_count - num_tested);
      GetRNGstate();
      a_select = runif(0.0, 1.0);
      PutRNGstate();
      if (a_select < select_prob) {
        GetRNGstate();
        a_select = runif(0.0, 1.0);
        PutRNGstate();
        if (a_select < 0.5) {
          plus_minus = -1;
        } else {
          plus_minus = 1;
        }
        GetRNGstate();
        increment = runif((0.0 - new_pop[i * num_species + j]), (1.0 - new_pop[i * num_species + j]));
        PutRNGstate();
        new_pop[i * num_species + j] += damper * increment * plus_minus;
        num_changed += 1;
      }
      num_tested += 1;
    }
    if ((i < num_species) && (j < (num_species - 1))) {
      j += 1;
    } else {
      if (j == (num_species - 1)) {
        j = 0;
        i += 1;
      }
    }
  }


  // return values
  for (i = 0; i < (num_species * num_species); i++) {
    if (new_pop[i] < 0) {
      new_pop[i] = 0.05;
    }
    if (new_pop[i] > 1) {
      new_pop[i] = 1;
    }
    pops[i] = new_pop[i];
  }

}