#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> // Includes `M_PI`
#include "timing.h"

#define k_B         (1.38065e-23)   // Boltzmann constant in [J/K]
#define stickcoeff  1               // Sticky coefficient in [0,1]; defines how likely particle to stick on grain if collides
#define T_gas       15              // Temperature of the gas in [K]
#define T_grain     15              // Temperature of the grain in [K]
#define mass_H      (1.6738234e-27) // Mass of hydrogen atom in [kg]
#define E_bind_H    (5.14982e-21)   // Binding energy of hydrogen atom to the grain surface in [Joules]; translated from 373K;
#define E_surf_H    (3.957e-21)     // Diffusion barrier energy of hydrogen atom on the grain surface;
// olivine surface has diffusion energy of 3.957e-21 Joules; amorphous carbon has diffusion energy of 7.0495771e-21 Joules
// v_bind = 1.116661e+12, v_surf = 0.978832e+12
#define Ns          (2e18)          // The surface density of sites, or cells, on the grain surface, 2e18m^-2 for olivine
#define p_no_event  (1-(1e-3))          // the probability that a hydrogen atom neither desorbs nor diffuses in a given time step
#define numdens     (1e18)    // Number density of a hydrogen particle in [m^-3]; just made up numbers for now

int n = 0;
int NUM_STEPS = 0;
int N = 0;
int TIMING_ONLY = 0;

// global count for how much H2 has been generated
int H2_count = 0;
int H_count = 0;
int H_desorption_count = 0;
double r_grain;  // Grain radius in [m]; from ((# of sites) = pi ( (grain radius*2)^2) * (surf. dens. of sites))
double ads_rate; // Calculate adsorption rate
double dt; // Typical time step 0.1 picoseconds (0.1e-12)

void initialize_occupancy(int occupancy[n][n]) {
  int i, j;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      occupancy[i][j] = 0;
}

double calc_rate(double E_process)
{
  return sqrt(2*Ns*E_process/(M_PI*M_PI*mass_H))*exp(-E_process/(k_B*T_grain));
}

void calc_probabilities(double buf[2])
{
// function to find movement probabilities and update matrix
  double R_diff = calc_rate(E_surf_H);
  double R_des = calc_rate(E_bind_H);

  // Probability of desorption
  buf[0] = (double) (1 - exp(-(R_diff + R_des)*dt))* R_des/ (R_diff + R_des);

  // probability of diffusion in any direction
  buf[1] = (double) (1 - exp(-(R_diff + R_des)*dt)) * R_diff/ ((R_diff + R_des));
}

// this function randomly samples from hop probabilities and updates the global occupancy matrix
void hop(int occupancy[n][n]) {
  int i,j;
  // indices: north, south, east, west, northwest, northeast, southwest, southeast;

  int directions[8][2] = {{0,1},{0,-1},{1,0},{-1,0},{-1,1},{1,1},{-1,-1},{1,-1}};

  // indices: 0:desportion, 1:diffusion
  double probs[2] = {0.0,0.0};

  calc_probabilities(probs);

  // prob_helper is cumulative probability
  double prob_helper[2] = {0.0,0.0};

  prob_helper[0] = probs[0];
  prob_helper[1] = probs[0] + probs[1];

  for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
      // cell is occupied
      if (occupancy[i][j] == 1) {
        double sample = random() / (double) RAND_MAX;

        // desportion
        if (sample < prob_helper[0]) {
          H_desorption_count++;
          occupancy[i][j] = 0;
        }

        // diffusion
        else if (sample > prob_helper[0] && sample < prob_helper[1]) {
          int random_index = random() % 8;
          int k = i + directions[random_index][0];
          int z = j + directions[random_index][1];
          // Deal with boundaries of the grain
          if (k < 0) { // Roll over the left
            k = n - 1;
          }
          if (k == n) { // Roll over the rightcd
            k = 0;
          }
          if (z < 0) { // Roll over the bottom
            z = n - 1;
          }
          if (z == n) { // Roll over the top
            z = 0;
          }
          if (occupancy[k][z] == 0) {
            occupancy[i][j] = 0;
            occupancy[k][z] = 1;
          }
          else {
            H2_count++;
            occupancy[k][z] = 0;
          }
        }
        // remaining in the same place is the complement
      }
    }
  }
}

double calc_timestep()
{
  double R_diff = calc_rate(E_surf_H);
  double R_des = calc_rate(E_bind_H);

  //p_no_event = 1 - e^(-(R_diff + R_des)*dt)
  return fmin(-log(p_no_event)/(R_diff + R_des), -log(p_no_event)/(ads_rate));
}

double calc_velocity()
{
  return sqrt(8*k_B*T_gas/(M_PI*mass_H));
}




double calc_next_adsorption_time (double last_adsorption_time) // Return next timestep upon which to add particle
{

  // Calculate and return next time this event will occur
  double randn = random() / (double) RAND_MAX; // Random number from 0 to 1
  return (-1 * log(randn) / ads_rate) + last_adsorption_time; // Random next timestep at which event will occur
}



void run_sim (int occupancy[n][n]) // Run full simulation for given number of timesteps
{
  double next_adsorption_time = 0.0; // Next timestep at which a hydrogen particle would be added

  // Set up loop through timesteps
  for (int ti=0; ti<NUM_STEPS; ti++) {
    if (next_adsorption_time <= ti * dt) { // If next adsorption timestep for this particle <= current time step
      // Choose a random place to deposit this particle
      int randx = (int) round((random() / (double) RAND_MAX) * n); // Random x-location
      int randy = (int) round((random() / (double) RAND_MAX) * n); // Random y-location
      H_count++; // Increment total count of added H atoms
      if (occupancy[randx][randy] == 1) {
        H2_count++;
        occupancy[randx][randy] = 0;
      }
      else {
        occupancy[randx][randy] = 1; // Add new hydrogen to occupancy matrix
      }

      // Update value in next_adsorption_time with the next timestep upon which to add this type of particle
      next_adsorption_time = calc_next_adsorption_time(next_adsorption_time);
    }

    // Loop through existing particles; determine movement and any reactions
    hop(occupancy);

  }
}


int main (int argc, char *argv[])
{
  if (argc == 4) {
    n = atoi(argv[1]);
    NUM_STEPS = atoi(argv[2]);
    TIMING_ONLY = atoi(argv[3]);
    N = n * n;
  }
  else {
    printf("arguments: n, NUM_STEPS, TIMING_ONLY\n");
    exit(-1);
  }

  timing_t tstart, tend;
  // start timing
  get_time(&tstart);

  // generate random seed
  srandom(time(NULL));

  int i, j, k;

  int occupancy[n][n];

  //fill global variables
  r_grain = (sqrt(N/1.0/(M_PI*Ns))/2.0);  // Grain radius in [m]; from ((# of sites) = pi ( (grain radius*2)^2) * (surf. dens. of sites))
  ads_rate = stickcoeff * calc_velocity() * numdens * M_PI * (r_grain * r_grain); // Calculate adsorption rate
  dt = calc_timestep(); // Typical time step 0.1 picoseconds (0.1e-12)

  initialize_occupancy(occupancy);
  run_sim(occupancy);

  // end timing
  get_time(&tend);
  if (TIMING_ONLY == 0) {
    printf("*****************************************************\n");
    // Print out result matrix...
    printf("Result Matrix:\n");
    for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
        printf("%d   ", occupancy[i][j]);
      printf("\n");
    }
    printf("******************************************************\n");
    printf ("Done.\n");
    printf("Elapsed time: %f s\n", (double)timespec_diff(tstart, tend));
    printf("Time simulated: %e s\n", dt*NUM_STEPS);
    printf("H desorption: %d\n", H_desorption_count);
    printf("H added: %d\n", H_count);
    printf("H2 generated: %d\n", H2_count);
    printf("Recombination efficiency (full, not past steady state): %f\n", (2*H2_count/1.0/H_count));
  }
  else if (TIMING_ONLY == 1) {
    printf("%f, ", (double)timespec_diff(tstart, tend));
  }
}
