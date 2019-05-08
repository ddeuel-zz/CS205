#include <math.h> // Includes `M_PI`
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "timing.h"

#define k_B         (1.38065e-23)   // Boltzmann constant in [J/K]
#define stickcoeff  1               // Sticky coefficient in [0,1]; defines how likely particle to stick on grain if collides
#define T_gas       50             // Temperature of the gas in [K]
#define T_grain     10             // Temperature of the grain in [K]
#define mass_H      (1.6738234e-27) // Mass of hydrogen atom in [kg]
#define E_bind_H    (5.14982e-21)   // Binding energy of hydrogen atom to the grain surface in [Joules]; translated from 373K;
#define E_surf_H    (3.957e-21)     // Diffusion barrier energy of hydrogen atom on the grain surface;
// olivine surface has diffusion energy of 3.957e-21 Joules; amorphous carbon has diffusion energy of 7.0495771e-21 Joules
// v_bind = 1.116661e+12, v_surf = 0.978832e+12
#define Ns          (2e18)          // The surface density of sites, or cells, on the grain surface, 2e18m^-2 for olivine
#define p_no_event  (1-(1e-3))          // the probability that a hydrogen atom neither desorbs nor diffuses in a given time step
#define numdens     (1e18)    // Number density of a hydrogen particle in [m^-3]; just made up numbers for now

// MPI_Datatype QUAD;
// #define QUADRANT(Q,groupy,groupx,subwidth) (Q[groupy * subwidth] + (groupx * subwidth))

int n = 0;
int nloc = 0;
int NUM_STEPS = 0;
int N = 0;
int csv_output = 0;

// global count for how much H2 has been generated
int H2_count = 0;
int H_count = 0;
int H_desorption_count = 0;
int H2_count_partial = 0;
int H_count_partial = 0;
int H_desorption_count_partial = 0;
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

void ghost_exchange(int occ[nloc][nloc], int colindingroup, int rowindingroup, int numcolingroup, int numrowingroup, int numgroupcolingrid, int numgrouprowingrid, int rank)
{
	// Note: Assumes rank is read left-right, up-down
	// Variables: colindingroup = index of column within this submatrix; rowindingroup = index of row within this submatrix
	//			: numcolingroup = number of columns within this submatrix; numrowingroup = number of rows within this submatrix; THESE DO NOT INCLUDE GHOST REGIONS
	//			: numgroupcolingrid = number of columns containing submatrices across entire grid (i.e., how many submatrices fit into grid along the horizontal axis); numgrouprowingrid = number of rows containing submatrices across entire grid (i.e., how many submatrices fit into grid along the vertical axis)

	if ((colindingroup == (numcolingroup-2)) || (colindingroup == (numcolingroup-1))) // Second-to-last OR last spot in row within this group
	{
		MPI_Barrier(MPI_COMM_WORLD); // Stop and wait at second-to-last OR last place in current row; shouldn't continue until all threads reach this spot
		
		// Send information FROM RIGHT GHOST REGION to the right
		if ((rank+1) % numgroupcolingrid == 0) // If at far-right of grid
		{
			MPI_Send(&occ[rowindingroup][nloc-1], 1, MPI_INT, (rank+1-numgroupcolingrid), 0, MPI_COMM_WORLD); // Info must roll over the edge of the grain and be sent to the far-left
		}
		else // If not at far-right of grid
		{
			MPI_Send(&occ[rowindingroup][nloc-1], 1, MPI_INT, (rank+1), 0, MPI_COMM_WORLD); // Info must be sent to the group to the right
		}
		occ[rowindingroup][nloc-1] = 0; // Clear after sending
		
		// Send information FROM LEFT GHOST REGION to the left
		if (rank % numgroupcolingrid == 0) // If at far-left of grid
		{
			MPI_Send(&occ[rowindingroup][0], 1, MPI_INT, (rank+numgroupcolingrid-1), 0, MPI_COMM_WORLD); // Info must roll over the edge of the grain and be sent to the far-right
		}
		else // If not at far-left of grid
		{
			MPI_Send(&occ[rowindingroup][0], 1, MPI_INT, (rank-1), 0, MPI_COMM_WORLD); // Info must be sent to the left
		}
		occ[rowindingroup][0] = 0; // Clear after sending
		
		// Receive information from the right INTO ACTUAL REGION OF DATA (NOT GHOST REGION)
		if ((rank+1) % numgroupcolingrid == 0) // If at far-right of grid
		{
			MPI_Recv(&occ[rowindingroup][numcolingroup], 1, MPI_INT, (rank+1-numgroupcolingrid), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must roll over from the edge of the grain and be received from far-left
		}
		else
		{
			MPI_Recv(&occ[rowindingroup][numcolingroup], 1, MPI_INT, (rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must be received from the right
		}
		
		// Receive information from the left INTO ACTUAL REGION OF DATA (NOT GHOST REGION)
		if (rank % numgroupcolingrid == 0) // If at far-left of grid
		{
			MPI_Recv(&occ[rowindingroup][1], 1, MPI_INT, (rank+numgroupcolingrid-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must roll over from edge of the grain and be received from far-right
		}
		else
		{
			MPI_Recv(&occ[rowindingroup][1], 1, MPI_INT, (rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must be received from the group to the left
		}
		
		
		// If not along bottom row of group, then return
		if (rowindingroup != numrowingroup)
			return;
		
		
		// Otherwise, if along bottom row of group, then exchange information
		
		// Send information to the top FROM TOP GHOST REGION
		if (rank < numgroupcolingrid) // If at far-top of the grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Send(&occ[0][jtemp], (nloc-2), MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Info must roll over edge of grain to be sent to far-bottom			
				occ[0][jtemp] = 0; // Clear after sending
			}
		}
		else // If not at far-top of grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Send(&occ[0][jtemp], (nloc-2), MPI_INT, (rank-numgroupcolingrid), 0, MPI_COMM_WORLD); // Info must be sent to group above
				occ[0][jtemp] = 0; // Clear after sending
			}
		}
		
		// Send information to the bottom FROM BOTTOM GHOST REGION
		if (rank >= (numgroupcolingrid*(numgrouprowingrid-1))) // If at far-bottom of grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Send(&occ[nloc-1][jtemp], (nloc-2), MPI_INT, (rank - ((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Info must roll over edge of grain to be sent to far-top		
				occ[nloc-1][jtemp] = 0; // Clear after sending
			}	
		}
		else // If not at far-bottom of grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Send(&occ[nloc-1][jtemp], (nloc-2), MPI_INT, (rank+numgroupcolingrid), 0, MPI_COMM_WORLD); // Info must be sent to bottom
				occ[nloc-1][jtemp] = 0; // Clear after sending
			}
		}
		
		// Receive information from the top INTO ACTUAL REGION OF DATA (NOT GHOST REGION)
		if (rank < numgroupcolingrid) // If at far-top of the grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Recv(&occ[1][jtemp], (nloc-2), MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must roll over edge of grain and be received from the group at the far-top
			}
		}
		else // If not at far-top of the grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Recv(&occ[1][jtemp], (nloc-2), MPI_INT, (rank-numgroupcolingrid), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must be received from the group at the top
			}
		}
		
		// Receive information from the bottom INTO ACTUAL REGION OF DATA (NOT GHOST REGION)
		if (rank >= (numgroupcolingrid*(numgrouprowingrid-1))) // If at far-bottom of grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Recv(&occ[numrowingroup][jtemp], (nloc-2), MPI_INT, (rank-(numgroupcolingrid*(numgrouprowingrid-1))), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must roll over edge of grain be received from the far-top
			}
		}
		else // If not at far-bottom of grain
		{
			int jtemp;
			for (jtemp=1; jtemp<(nloc-1); jtemp++)
			{
				MPI_Recv(&occ[numrowingroup][jtemp], (nloc-2), MPI_INT, (rank+numgroupcolingrid), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Info must be received from the bottom
			}
		}
		
		// Deal with obnoxious ghost-cell corner cases
		if (rank >= (numgroupcolingrid*(numgrouprowingrid-1))) // If at far-bottom of grain
		{
			// Send to actual data at top-left and top-right
			MPI_Send(&occ[nloc-1][0], 1, MPI_INT, (rank - ((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Deals with pesky bottom-left ghost cell; info must roll over edge of grain to be sent to top-right cell			
			occ[nloc-1][0] = 0; // Clear after sending
			MPI_Send(&occ[nloc-1][nloc-1], 1, MPI_INT, (rank - ((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Deals with pesky bottom-right ghost cell; info must roll over edge of grain to be sent to top-left cell			
			occ[nloc-1][nloc-1] = 0; // Clear after sending
			
			// Recieve from ghost cells at top-left and top-right
			MPI_Recv(&occ[numrowingroup][nloc-2], 1, MPI_INT, (rank-(numgroupcolingrid*(numgrouprowingrid-1))), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Deals with pesky bottom-right ghost cell; info must roll over edge of grain to be received from top-left cell
			MPI_Recv(&occ[numrowingroup][1], 1, MPI_INT, (rank-(numgroupcolingrid*(numgrouprowingrid-1))), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Deals with pesky bottom-left ghost cell; info must roll over edge of grain to be received from top-right cell
		}
		
		if (rank < numgroupcolingrid) // If at far-top of the grain
		{
			// Send to actual data at bottom-left and bottom-right
			MPI_Send(&occ[0][0], 1, MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Deals with pesky top-left ghost cell; info must roll over edge of grain to be sent to bottom-right cell
			occ[0][0] = 0; // Clear after sending
			MPI_Send(&occ[0][nloc-1], 1, MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD); // Deals with pesky top-right ghost cell; info must roll over edge of grain to be sent to bottom-left cell
			occ[0][nloc-1] = 0; // Clear after sending
			
			// Receive from ghost cells at bottom-left and bottom-right
			MPI_Recv(&occ[1][nloc-2], 1, MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Deals with pesky top-right ghost cell; info must roll over edge of grain to be received from bottom-left cell
			MPI_Recv(&occ[1][1], 1, MPI_INT, (rank+((numgrouprowingrid-1)*numgroupcolingrid)), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Deals with pesky top-left ghost cell; info must roll over edge of grain to be received from bottom-right cell
		}
	}
}


// this function randomly samples from hop probabilities and updates the local occ matrix
void hop(int occ[nloc][nloc], int nper, int rank, int size) {
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
 
  for (i=1; i<nloc-1; i++) { // NOTE: Padded by +1 so as not to include ghost region
   for (j=1; j<nloc-1; j++) { // NOTE: Padded by +1 so as not to include ghost region
      // cell is occupied
      if (occ[i][j] == 1) {
        double sample = random() / (double) RAND_MAX;

        // desorption
        if (sample < prob_helper[0]) {
          H_desorption_count_partial++;
          occ[i][j] = 0;
        }

        // diffusion
        else if (sample > prob_helper[0] && sample < prob_helper[1]) {
          int random_index = random() % 8;
          int k = i + directions[random_index][0];
          int z = j + directions[random_index][1];
          // Deal with boundaries of the grain
	      // NOTE: NONE OF THIS SHOULD EVER BE CALLED; TAKEN CARE OF IN GHOST_EXCHANGE FUNCTION
          if (k < 0) { // Roll over the left
            // k = n - 1;
			printf("ERROR: SOMETHING WENT WRONG WITH PERUSING SUBMATRIX. MAKE SURE THAT GHOST REGIONS ARE NEVER REACHED.");
			exit(-1);
          }
          if (k == nloc) { // Roll over the rightcd
            // k = 0;
			printf("ERROR: SOMETHING WENT WRONG WITH PERUSING SUBMATRIX. MAKE SURE THAT GHOST REGIONS ARE NEVER REACHED.");
			exit(-1);
          }
          if (z < 0) { // Roll over the bottom
            // z = nloc - 1;
			printf("ERROR: SOMETHING WENT WRONG WITH PERUSING SUBMATRIX. MAKE SURE THAT GHOST REGIONS ARE NEVER REACHED.");
			exit(-1);
          }
          if (z == nloc) { // Roll over the top
            // z = 0;
			printf("ERROR: SOMETHING WENT WRONG WITH PERUSING SUBMATRIX. MAKE SURE THAT GHOST REGIONS ARE NEVER REACHED.");
			exit(-1);
          }
		  
		  // Exchange info gathered in ghost regions of submatrix
		   ghost_exchange(occ, i, j, nper, nper, size, size, rank);
		   
		  // Place particle; check for reaction 
          if (occ[k][z] == 0) {
            occ[i][j] = 0;
            occ[k][z] = 1;
          }
          else {
            H2_count_partial++;
            occ[k][z] = 0;
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



void run_sim (int occ[nloc][nloc], int nper, int rank, int size) // Run full simulation for given number of timesteps
{
  double next_adsorption_time = 0.0; // Next timestep at which a hydrogen particle would be added

  // Set up loop through timesteps
  for (int ti=0; ti<NUM_STEPS; ti++) {
    if (next_adsorption_time <= ti * dt) { // If next adsorption timestep for this particle <= current time step


      // Choose a random place to deposit this particle
	int randrank = (int) round((random() / (double) RAND_MAX) * size); // Random submatrix within which to place particle
	if (randrank >= size)
		exit(-1);
	if (rank == randrank) // If this rank's section contains particle
	{
         int randx = (int) round((random() / (double) RAND_MAX) * nper); // Random x-location
         int randy = (int) round((random() / (double) RAND_MAX) * nper); // Random y-location
	      H_count_partial++; // Increment total count of added H atoms
	      if (occ[randx+1][randy+1] == 1) {
		H2_count_partial++; // NOTE: Padded by +1 so as not to include ghost region
		occ[randx+1][randy+1] = 0; // NOTE: Padded by +1 so as not to include ghost region
	      }
	      else {
		occ[randx+1][randy+1] = 1; // Add new hydrogen to occ matrix; NOTE: Padded by +1 so as not to include ghost region
	      }

	      // Update value in next_adsorption_time with the next timestep upon which to add this type of particle
	      next_adsorption_time = calc_next_adsorption_time(next_adsorption_time);
	    }
	}

    // Loop through existing particles; determine movement and any reactions
    hop(occ, nper, rank, size);
  }
}


int main (int argc, char *argv[])
{
  if (argc == 4) {
    n = atoi(argv[1]);
    NUM_STEPS = atoi(argv[2]);
    csv_output = atoi(argv[3]);
    N = n * n;
  }
  else {
    printf("arguments: n, NUM_STEPS, csv_output\n");
    exit(-1);
  }
  

  timing_t tstart, tend;
  double tstartMPI, tendMPIsim, tendMPIall;
  // start timing
  get_time(&tstart);

  // generate random seed
  srandom(10);

  int i, j, k;

  int occupancy[n][n];

  //fill global variables
  r_grain = (sqrt(N/1.0/(M_PI*Ns))/2.0);  // Grain radius in [m]; from ((# of sites) = pi ( (grain radius*2)^2) * (surf. dens. of sites))
  ads_rate = stickcoeff * calc_velocity() * numdens * M_PI * (r_grain * r_grain); // Calculate adsorption rate
  dt = calc_timestep(); // Typical time step 0.1 picoseconds (0.1e-12)

  initialize_occupancy(occupancy); // Initialize occupancy

  
  // Initialize MPI; get rank and size
  int rank, size, nper, ioffset, joffset;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  tstartMPI = MPI_Wtime();
  

  // Print diagnostic message
  if (rank == 0)
	printf("Processes: %d\n", size);
	
  // Split up occupancy between processes
  nper = (int) sqrt(n*n / size); // Allocate size of submatrix per processor
  if (n % nper != 0) // Raise error if unequal distribution of matrix between processors
  {
    printf("ERROR: Please choose grain size so that it can be split evenly, in both length and width, among desired number of processors.");
	exit(-1);
  }
  ioffset = (rank / (n / nper))*nper; // Starting point along y-axis in matrix for current processor
  joffset = (rank - ((rank / (n / nper)) * (n / nper)))*nper; // Starting point along x-axis in matrix for current processor
  nloc = nper + 2; // Length=width for each submatrix; padded by 1 to left, 1 to right, 1 to top, 1 to bottom for ghost_exchange functionality
    

  int itemp, jtemp;
  int occpiece[nloc][nloc];
  for (itemp=0; itemp < nloc; itemp++)
  {
    for (jtemp=0; jtemp < nloc; jtemp++)
    {
        occpiece[itemp][jtemp] = 0;
    }
  }


  run_sim(occpiece, nper, rank, size); // Run simulation for this submatrix
    tendMPIsim = MPI_Wtime();

  if (csv_output == 0) {  
  {
	int pind, temprow, tempcol;
	// Store results from all processes into the occupancy matrix
	if (rank == 0) // Fill out with current rank=0 as a base
	{
	    // Fill data from rank=0 submatrix
	    for (i=1; i<(nloc-1); i++)
	    {
	       for (j=1; j<(nloc-1); j++)
	       {
			occupancy[i-1][j-1] = occpiece[i][j];
		}
	    }
	    // Receive and then store data from every other matrix
	    for (pind=1; pind<size; pind++)
	    {
		for (i=1; i<(nloc-1); i++)
		{
			MPI_Recv(&occpiece[i], nper, MPI_INT, pind, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive info from process pind along this row
		}
		temprow = pind / (n / nper); // This section's row
		tempcol = pind - (temprow*(n / nper)); // This section's column
		// Fill data from the extracted submatrix from pind
		for (i=1; i<(nloc-1); i++)
		{
		       for (j=1; j<(nloc-1); j++)
		       {
				occupancy[i-1+(temprow*nper)][j-1+(tempcol*nper)] = occpiece[i][j];
			}
		 }
	    }
	}
	else // Send info from this rank (for rank != 0)
	{
		for (i=1; i<(nloc-1); i++)
		{
		MPI_Send(&occpiece[i], nper, MPI_INT, 0, 0, MPI_COMM_WORLD); // Send info to rank 0
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&H_count_partial, &H_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&H_desorption_count_partial, &H_desorption_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&H2_count_partial, &H2_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    // Below prints results, if rank==0
   /* if (rank == 0)
    {
	    printf("*****************************************************\n");
	    printf("Result Matrix:\n"); // Receive info from processors to print
	    for (i=0; i<n; i++)
	    {
		for (j=0; j<n; j++)
		{
		  printf("%d  ", occupancy[i][j]);
		}
    		printf("\n");
	    }
	printf("*****************************************************\n"); 
    } */
  }
  // end timing
  get_time(&tend);
  tendMPIall = MPI_Wtime();
   if (rank == 0)
    {
    printf ("Done.\n");
    printf("Elapsed time: %f s\n", (double)timespec_diff(tstart, tend));
    printf("Elapsed MPI sim time: %f s\n", (double)(tendMPIsim - tstartMPI));
    printf("Elapsed MPI all time: %f s\n", (double)(tendMPIall - tstartMPI));
    printf("Time simulated: %e s\n", dt*NUM_STEPS);
    printf("H desorption: %d\n", H_desorption_count);
    printf("H added: %d\n", H_count);
    printf("H2 generated: %d\n", H2_count);
    printf("Recombination efficiency (full, not past steady state): %f\n", (2*H2_count/1.0/H_count));
    }
   }
  else if (csv_output == 1) {
    if (rank == 0)
       printf("%Lf, %d, %d", (double)timespec_diff(tstart, tend), n, NUM_STEPS);
  }
    MPI_Finalize();
}
