/*
 *  Dominance calculation for food web genetic algorithm
 *  Author: Jian Yen
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>

void GA_mutate_select(double *pops, int *dom, int *n_pop, int *n_keep, int *n_species,
                      int *k_val, double *k_prob)
{
  // Initialise
  int i, j, k;
  int num_pop = *n_pop;
  int num_keep = *n_keep;
  int num_species = *n_species;
  int k_vals = *k_val;
  double k_prob_set = *k_prob;
  double k_probs[k_vals];
  double *new_pop;
  new_pop = malloc (num_species * num_species * sizeof(double));

  int zero_count = 0, ones_count = 0;
  int change_count = 0;
  int change_range;
  double damper = 1;
  int num_changed = 0;
  int num_tested = 0;
  int i_counter = num_keep;
  int plus_minus;
  int num_to_change;
  
  int k_cur = 0;
  int choose_pop[k_vals];
  double select_prob;
  double a_select;

  int swapped = 1;
  int temp;
  int pop_select = 0;

  int new_pop_test = 0;
  int new_pop_clash = 1;
  
  // Store pops and dom in appropriate arrays
  int dom_vals[num_pop];
  double *pops_all;
  pops_all = malloc (num_pop * num_species * num_species * sizeof(double));

  for (i = 0; i < num_pop; i++) {
    dom_vals[i] = dom[i];
    for (j = 0; j < num_species; j++) {
      for (k = 0; k < num_species; k++) {
        pops_all[i * num_species * num_species + j * num_species + k] = pops[i * num_species * num_species + j * num_species + k];
      }
    }
  }

  k_probs[0] = k_prob_set;
  for (i = 1; i < k_vals; i++) {
    k_probs[i] = k_probs[i - 1] + k_prob_set * pow((1 - k_prob_set), i);
  }
  for (i = 0; i < k_vals; i++) {
    k_probs[i] /= k_probs[(k_vals - 1)];
  }  

  // Run main selection and mutation
  while (i_counter < num_pop) {

    // pull out k_vals populations at random
    for (i = 0; i < k_vals; i++) {
      choose_pop[i] = 0;
    }
    i = 0;
    k_cur = 0;
    while (k_cur < k_vals) {
      select_prob = (double)(k_vals - k_cur) / (double)(num_keep - i);
      GetRNGstate();
      a_select = runif(0.0, 1.0);
      PutRNGstate();
      if (a_select < select_prob) {
        choose_pop[k_cur] = i;
        k_cur += 1;
      }
      i += 1;
    }

    // set up a bubble sort/rank to put choose_pop in order of increasing dom_vals
    while (swapped) {
      swapped = 0;
      for (i = 0; i < k_vals; i++) {
        for (j = 0; j < (k_vals - 1); j++) {
          if (dom_vals[choose_pop[i]] > dom_vals[choose_pop[i + 1]]) {
            temp = choose_pop[i];
            choose_pop[i] = choose_pop[i + 1];
            choose_pop[i + 1] = temp;
            swapped = 1;
          }
        }
      }
    }
   
    // sample from choose_pop with probability k_probs
    GetRNGstate();
    a_select = runif(0.0, 1.0);
    PutRNGstate();
    pop_select = 0;
    while (a_select > k_probs[pop_select]) {
      pop_select += 1;
    }
       
    // set new_pop based on choose_pop[i]
    zero_count = 0;
    ones_count = 0;
    for (j = 0; j < num_species; j++) {
      for (k = 0; k < num_species; k++) {
        new_pop[j * num_species + k] = pops_all[choose_pop[pop_select] * num_species * num_species + j * num_species + k];
        if (new_pop[j * num_species + k] == 0) {
          zero_count += 1;
        } else {
          ones_count += 1;
        }
      }
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
      if (pops_all[choose_pop[pop_select] * num_species * num_species + i * num_species + j] != 0.0) {
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
          new_pop[i * num_species + j] += damper * plus_minus;
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

    // set negative values to small positive and values > 1 to 1
    for (i = 0; i < num_species; i++) {
      for (j = 0; j < num_species; j++) {
        if (new_pop[i * num_species + j] < 0.0) {
          new_pop[i * num_species + j] = 0.05;
        }
        if (new_pop[i * num_species + j] > 1) {
          new_pop[i * num_species + j] = 1;
        }
      }
    }

    // check for duplicates
    i = 0;
    while ((i < i_counter) && (new_pop_clash == 1)) {
      new_pop_test = 0;
      for (j = 0; j < num_species; j++) {
        for (k = 0; k < num_species; k++) {
          if (new_pop[j * num_species + k] == pops_all[i * num_species * num_species + j * num_species + k]) {
            new_pop_test += 1;
          }
        }
      }
      if (new_pop_test == (num_species * num_species)) {
        new_pop_clash = 0;
      }
      i += 1;
    }

    if (new_pop_clash != 0) {
      for (j = 0; j < num_species; j++) {
        for (k = 0; k < num_species; k++) {
          pops_all[i_counter * num_species * num_species + j * num_species + k] = new_pop[j * num_species + k];
        }
      }
      i_counter += 1;
      damper *= 0.8;
    }
  }

  //Update and return values
  for (i = 0; i < num_pop; i++) {
    dom[i] = dom_vals[i];
    for (j = 0; j < num_species; j++) {
      for (k = 0; k < num_species; k++) {
        pops[i * num_species * num_species + j * num_species + k] = pops_all[i * num_species * num_species + j * num_species + k];
      }
    }
  }
  
  free(pops_all);
  
}