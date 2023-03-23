
#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <sstream>

using namespace std;

/* This struct contains parameters that are used by the model. They can be specified as commandline parameters.
 * Currently only randseed, foodinglucose, celldrift, and o2 are specified as commandline parameters
 * the rest keeps the default value.
 */

Par::Par()
{
  /* Construct parameter struct */
  datadir = "/data";                 // directory in which data is placed
  moviedir = "/movie";               // directory in which png-files are placed
  movie = true;                      // whether simulation produces .png files to make a movie
  timescale = 20.0;                  // every timestep is 1 hr/20 = 3 minutes
  cellfrac = 0.025;                  // OVERWRITTEN
  MIN_GROWTH = 0.0000001;            // Minimum growth to bother having fluxes at all. If solution gives a growth below this, just set all fluxes to 0
  DIFF_CUTOFF = 0.0000001;           // Don't bother diffusion substances with a concentration beneath this
  INITDENS = 0.3;                    // Initial cell density, fraction of empty squares that will get a bacterium if not using bacplace
  foodinglucose = 0;                 // OVERWRITTEN BY RUN_ALL-default amount of glucose that is given on every pulse umol per feeding, altered by command
  foodingalactose = 0;               // OVERWRITTEN galactose in umol per feeding
  foodinlactose = 0;                 // OVERWRITTEN lactose in umol per feeding
  foodinlactate = 0;                 // OVERWRITTEN lactate in umol per feeding
  foodinGOS = 0;                 // OVERWRITTEN lactate in umol per feeding
  foodinfl = 0;                 // OVERWRITTEN lactate in umol per feeding
  timefood = 3;                      // amount of time between pulses: 3hr. Can be put to some implausibly high number to only have food at time 0
  divisionsize = 2.0;                // Size at which to divide
  diffconst = 2;                     // OVERWRITTEN
  pdeath = 2.0;                      // OVERWRITTEN random death rate of cells (density dependent): per hour
  pdeath_base = 0.08;                // OVERWRITTEN basal random death rate: per hour
  ncol = 225;                        // grid size (length)
  nrow = 8;                          // grid size (width)
  metdrift = 20;                     // Metabolic drift (if set to 10, move every other step, if set to 20, move every step. If set to 40, move twice every step)
  celldrift = 0;                     // OVERWRITTEN BY RUN_ALL-whether cells also flow: 0 no flow, 1 flow as fast as metabolites
  cellinflux = 0;                    // probability that cell that flows out at the end enters at the beginning
  shuffle = false;                   // mixes all cells and metabolites every time-step
  presolve = false;                  // presolving of metabolic decreases computational efficiency
  maxcells = 1;                      // maximum number of cells in one grid point
  randseed = 1000;                   // OVERWRITTEN BY RUN_ALL-random seed used
  runid = 1;                         // Identifies otherwise identical runs
  max_metup = 10000000000;           // maximal uptake rate
  bacplacetype = 2;                  // 2;//OVERWRITTEN BY RUN_ALL-place all bacteria on initialisation (1) or only some (0)
  bacplacefile = "bacplace";         // Default text file for loading in bacteria
  immigration = 2.0001;              // OVERWRITTEN BY RUN_ALL First number is immigration setting, second immigration chance per eligible grid square per step
  deathsetting = 2;                  // OVERWRITTEN BY RUN_ALL How death works, based on presets or a flat chance 0<n<1
  placechance = 0.0001;              // as the uniform random number generator is a six digit number, lowest possible value >0 is 0.000001
  growthmod = 2;                     // OVERWRITTEN Modifier to increase or decrease grow speed by making the biomass reaction cheaper or more expensive
  bactmove = 0;                      // OVERWRITTEN Determines whether bacteria move even when already feeding
  fluxweight = 0.00025;              // OVERWRITTEN .00025;//flux weight for all cells, normally 0.00025. Total flux allowed, if all fw equal, is equal to 0.2/fw. So, fw is allowed to sum to 0.2.
  testname = "blank";                // Name of experiment
  initialox = 0.1;                   // OVERWRITTEN BY RUN_ALL
  oxin = 0;                          // Oxygen influx per step
  biomassreplace = true;             // OVERWRITTEN
  maxswaps = 100000;                 // How many kawasaki steps (simulated mixing) we should do each step at most. This number ensures we will do the maximum amount
  initialvolume = 1;                 // 0.03; //Initial volume (size) of each population //OVERWRITTEN
  baseintake = 1;                    // Can be used to decrease initial lactose amount
  intakeincrease = 0;                // Can be used to gradually increase lactose input
  frequentdata = false;              // Save data frequently
  knockoutlactate = false;           // Knock out lactate uptake
  knockoutxfp = false;               // Knock out XFP in bifids
  knockoutnonbiflactose = false;     // Knock out lactate uptake for non-bifids
  knockoutnonbiflactoselate = false; // Knock out lactose uptake for non-bifids after 5040 steps
  knockoutbutyrolactoselate = false; //Knock out lactose uptake by butyrogenics after 5040 steps
  knockoutenterooxygen = false;      // Knock out oxygen uptake for e. coli
  knockoutbiflactate = false;        // Knock out lactate uptake for bifids
  knockoutbiflactateup = false;        // Knock out lactate uptake for bifids
  knockoutecollactate = false;       // Knock out lactate uptake for e. coli
  knockoutbifandecollactate = false; // Knock out lactate uptake in both e. coli and bifids
  knockoutbutyro12ppd=false;
  knockoutbutyrolcts=false;
  knockoutbutyrolactate=false;
  knockout12ppd=false;
  alwayssaveenergy = false;          // OVERWRRITTEN always save all energies, even if plausible
  diminishoxygenrelease=false;
  

  // Define the metabolites whose influx and outflux will be saved in exchange.dat
  outputmetabs.insert("gal");
  outputmetabs.insert("for");
  outputmetabs.insert("ac");
  outputmetabs.insert("for");
  outputmetabs.insert("ppa");
  outputmetabs.insert("succ");
  outputmetabs.insert("lac_D");
  outputmetabs.insert("lac_L");
  outputmetabs.insert("lcts");
  outputmetabs.insert("co2");
  outputmetabs.insert("h2");
  outputmetabs.insert("etoh");
  outputmetabs.insert("h2o");
  outputmetabs.insert("o2");
  outputmetabs.insert("but");
  outputmetabs.insert("2FuLa");
  outputmetabs.insert("GOS3");
  outputmetabs.insert("GOS4");
  outputmetabs.insert("GOS5");
outputmetabs.insert("12ppd_S");
}

Par::~Par()
{
  /* Destruct parameter struct */
}

Par par;
