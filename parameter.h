
#ifndef _PARAMETER_H_
#define _PARAMETER_H_
#include <vector>
#include <iostream>
#include <string>
#include <set>
using namespace std;

/* This struct contains parameters that are used by the model. They can be specified as commandline parameters.
 * Currently only randseed, foodinglucose, celldrift, and o2 are specified as commandline parameters
 * the rest keeps the default value.
 */

struct Par
{
  Par();                // Construct parameter struct
  ~Par();               // Destruct parameter struct
  std::string datadir;  // directory in which data is placed
  std::string moviedir; // directory in which png-files are placed
  bool movie;           // whether simulation produces .png files to make a movie
  double timescale;     // every timestep is 1 hr/20 = 3 minutes
  double cellfrac;      // 0.25-1 gr Dry weight/liter
  double MIN_GROWTH;
  double DIFF_CUTOFF;
  double INITDENS;      // inititial cell density: probability
  double foodinglucose; // amount of glucose that is given on every pulse mmol/liter
  double foodingalactose;
  double foodinlactose;
  double foodinlactate;
  double foodinGOS;
  double foodinfl;
  double oxin;
  int timefood;       // amount of time between pulses: 8hr
  double diffconst;   // diffusion constant of metabolites corresponds to 0.5*0.2cm*0.2cm/3600 =555 (micrometer)^2/s
  double pdeath;      // random death rate of cells (density dependent)
  double pdeath_base; // basal random death rate
  int ncol;           // grid size (x)
  int nrow;           // grid size (y)
  double metdrift;    // flow rate of metabolites through the gut
  double celldrift;   // whether cells also flow: 0 no flow, 1 flow as fast as metabolites
  double cellinflux;  // probability that cell that flows out at the end enters at the beginning
  bool shuffle;     // mixes all cells and metabolites every time-step
  bool presolve;    // presolving of metabolic decreases computational efficiency
  int maxcells;     // maximum number of cells in one grid point
  int randseed;     // random seed used
  int runid;        // to identify runs with the same name
  double max_metup; // maximal uptake rate
  int bacplacetype;
  std::string bacplacefile; //
  double immigration;       //
  double placechance;
  double growthmod;
  double deathsetting; // Death rate setting. 0>i>1 for a purely linear setting with death chance i, 2 for density-dependent
  int bactmove;
  double fluxweight;
  std::set<string> outputmetabs;
  std::string testname;
  double initialox;
  bool biomassreplace;
  int maxswaps;
  double initialvolume;
  double baseintake;     // Base intake of food on day 0
  double intakeincrease; // Extra intake based on step number
  bool frequentdata;     // Whether to store all data points, or only every 60th
  bool knockoutlactate;
  bool knockoutxfp;
  bool knockoutnonbiflactose;
  bool knockoutnonbiflactoselate;
  bool knockoutbutyrolactoselate;
  bool knockoutenterooxygen;
  bool knockoutbiflactate;
  bool knockoutbiflactateup;
  bool knockoutecollactate;
  bool knockoutbifandecollactate;
  bool knockoutbutyro12ppd;
  bool knockoutbutyrolcts;
  bool knockoutbutyrolactate;
  bool knockout12ppd;
  double divisionsize;
  bool alwayssaveenergy;
  bool diminishoxygenrelease;
};

#endif
