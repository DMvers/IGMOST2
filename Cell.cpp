#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Cell.h"
#include "parameter.h"

#define SIMPLEX 1 // interior gives sometimes numerical instability

/* This class implements the bacterium (metapopulation of bacteria) that occupies one position.
 * The class is constructed on the basis of an instance of class BactParam that contains the metabolic network for a specific species.
 * This class contains methods to use that metabolic network through FBA to grow and consume and produce metabolites.
 */

extern Par par;
extern int ancestorcount;
extern std::string simulationdir;

Cell::Cell()
{
  /* Method to construct an empty cell without a bacterium, so an empty position. */
  lp = glp_create_prob();
  params = (glp_smcp *)malloc(sizeof(glp_smcp));
  glp_init_smcp(params);

  species = "none";
  indexnr = 0;
  exist = false;
  bm_reaction = 0;

  nummetabolites = 0;
  numreactions = 0;
  numstoich = 0;

  exchange_dict.clear();
  exchange_list.clear();

  react_dict.clear();

  obj.clear();
  SparseS.clear();
  iRow.clear();
  jCol.clear();
  genome.clear();
  used.clear();

  lb.clear();
  ub.clear();
  flux.clear();
  fluxtot.clear();

  exist = false;
  hasmoved = false;
  active = false;
  exitflag = 0;
  ipos = 0;
  jpos = 0;
  timebirth = 0;
  time_alive = 0;
  growthrate = 0;
  volume = 0;
  ancestor = 0;
}

Cell::Cell(BactParam bp)
{
  /* Method to construct a Cell with a bacterium, the species is determined by the argument bp, which contains the metabolic network. */
  lp = glp_create_prob();
  params = (glp_smcp *)malloc(sizeof(glp_smcp));
  glp_init_smcp(params);

  species = bp.speciesName;
  indexnr = bp.indexnr;
  bm_reaction = bp.bm_reaction;

  nummetabolites = bp.nummetabolites;
  numreactions = bp.numreactions;
  numstoich = bp.numstoich;

  exchange_dict = bp.exchange_dict;
  exchange_list = bp.exchange_list; // exchange reactions

  maxflux = 0.2; // 0.75;//0.25;%as article bmc

  react_dict = bp.react_dict;

  obj = bp.objs;
  SparseS = bp.SparseS;
  iRow = bp.iRow;
  jCol = bp.jCol;
  genome = bp.genome;
  used = bp.used;

  lb = bp.lbs;
  ub = bp.ubs;
  flux = bp.fluxs;
  numreactions = bp.numreactions;
  numstoich = bp.numstoich;

  exchange_dict = bp.exchange_dict;
  exchange_list = bp.exchange_list; // exchange reactions

  maxflux = 0.2; // 0.75;//0.25;%as article bmc

  react_dict = bp.react_dict;

  obj = bp.objs;
  SparseS = bp.SparseS;
  iRow = bp.iRow;
  jCol = bp.jCol;
  genome = bp.genome;
  used = bp.used;

  lb = bp.lbs;
  ub = bp.ubs;
  flux = bp.fluxs;
  fluxtot.resize(numreactions);
  exist = false;
  hasmoved = false;
  active = false;
  exitflag = 0;
  ipos = 0;
  jpos = 0;
  timebirth = 0;
  time_alive = 0;

  growthrate = 0;

  volume = 0;
  ancestor = 0;
  fluxtot.resize(numreactions);
  exist = false;
  hasmoved = false;
  active = false;
  exitflag = 0;
  ipos = 0;
  jpos = 0;
  timebirth = 0;
  time_alive = 0;

  growthrate = 0;

  volume = 0;
  ancestor = 0;
}

Cell::Cell(const Cell &c)
{
  /* Method to construct a Cell on the basis of another instance of Cell. */
  lp = glp_create_prob();
  params = (glp_smcp *)malloc(sizeof(glp_smcp));
  glp_init_smcp(params);
  glp_copy_prob(lp, c.lp, 1);

  exchange_dict = c.exchange_dict;
  exchange_list = c.exchange_list; // exchange reactions
  react_dict = c.react_dict;

  genome = c.genome;
  used = c.used;
  obj = c.obj;
  SparseS = c.SparseS;
  iRow = c.iRow;
  jCol = c.jCol;
  lb = c.lb;
  ub = c.ub;
  lbubfood = c.lbubfood;
  flux = c.flux;
  fluxtot.resize(c.numreactions);

  genome.squeeze();
  used.squeeze();
  obj.squeeze();
  SparseS.squeeze();
  iRow.squeeze();
  jCol.squeeze();
  lb.squeeze();
  ub.squeeze();
  lbubfood.squeeze();
  flux.squeeze();
  fluxtot.squeeze();

  species = c.species;
  indexnr = c.indexnr;
  bm_reaction = c.bm_reaction;
  numreactions = c.numreactions;
  nummetabolites = c.nummetabolites;
  numstoich = c.numstoich;
  maxflux = c.maxflux;

  exchange_dict = c.exchange_dict;
  exchange_list = c.exchange_list; // exchange reactions
  react_dict = c.react_dict;

  genome = c.genome;
  used = c.used;
  obj = c.obj;
  SparseS = c.SparseS;
  iRow = c.iRow;
  jCol = c.jCol;
  lb = c.lb;
  ub = c.ub;
  lbubfood = c.lbubfood;
  flux = c.flux;
  fluxtot.resize(c.numreactions);

  genome.squeeze();
  used.squeeze();
  obj.squeeze();
  SparseS.squeeze();
  iRow.squeeze();
  jCol.squeeze();
  lb.squeeze();
  ub.squeeze();
  lbubfood.squeeze();
  flux.squeeze();
  fluxtot.squeeze();

  species = c.species;
  indexnr = c.indexnr;
  bm_reaction = c.bm_reaction;
  numreactions = c.numreactions;
  nummetabolites = c.nummetabolites;
  numstoich = c.numstoich;
  maxflux = c.maxflux;
  growthrate = c.growthrate;

  volume = c.volume;
  ancestor = c.ancestor;

  exist = c.exist;
  hasmoved = c.hasmoved;
  active = c.active;
  exitflag = c.exitflag;
  ipos = c.ipos;
  jpos = c.jpos;
  timebirth = c.timebirth;
  time_alive = c.time_alive;
  numancestors = c.numancestors;

  volume = c.volume;
  ancestor = c.ancestor;

  exist = c.exist;
  hasmoved = c.hasmoved;
  active = c.active;
  exitflag = c.exitflag;
  ipos = c.ipos;
  jpos = c.jpos;
  timebirth = c.timebirth;
  time_alive = c.time_alive;
  numancestors = c.numancestors;
}

Cell::~Cell()
{
  /* Destructor of Cell class. */
  glp_delete_prob(lp);
  free(params);
}

Cell &Cell::operator=(const Cell &c)
{
  /* Method to create a reference to the Cell class. */
  if (this == &c)
  {
    return *this;
  }
  exchange_dict = c.exchange_dict;
  exchange_list = c.exchange_list; // exchange reactions
  react_dict = c.react_dict;

  genome = c.genome;
  used = c.used;
  obj = c.obj;
  SparseS = c.SparseS;
  iRow = c.iRow;
  jCol = c.jCol;
  lb = c.lb;
  ub = c.ub;
  lbubfood = c.lbubfood;
  flux = c.flux;
  fluxtot.resize(c.numreactions);

  genome.squeeze();
  used.squeeze();
  obj.squeeze();
  SparseS.squeeze();
  iRow.squeeze();
  jCol.squeeze();
  lb.squeeze();
  ub.squeeze();
  lbubfood.squeeze();
  flux.squeeze();
  fluxtot.squeeze();

  species = c.species;
  indexnr = c.indexnr;
  bm_reaction = c.bm_reaction;
  numreactions = c.numreactions;
  nummetabolites = c.nummetabolites;
  numstoich = c.numstoich;
  maxflux = c.maxflux;

  growthrate = c.growthrate;

  volume = c.volume;
  ancestor = c.ancestor;

  exist = c.exist;
  hasmoved = c.hasmoved;
  active = c.active;
  exitflag = c.exitflag;
  ipos = c.ipos;
  jpos = c.jpos;
  timebirth = c.timebirth;
  time_alive = c.time_alive;
  numancestors = c.numancestors;

  return *this;
}

bool Cell::CanGrow() const
{
  // Determines whether FBA is performed, to speed-up simulations
  return true; // ignore, always perform FBA
}

void Cell::Do_Metabolism(QVector<BactParam> bps) // Ub of cell already set by Grid::ReadEnv()
{
  /* Cells grow on the metabolites available at their position.
   * This method calls the method for FBA, sees if a threshold for minimal growth is achieved and calculates total fluxes.
   * Mostly here for historic reasons, can be expanded to include other forms of metabolism (such as seperate polysaccharide metabolism),
   * could perhaps be joined with do_fba() otherwise.
   */
  do_fba(bps); // perform the fba
  if (growthrate < par.MIN_GROWTH || exitflag == 10)
  {
    growthrate = 0;
    for (int k = 0; k < numreactions; k++)
    {
      flux[k] = 0;
    }
  }

  if (exitflag != 0 && exitflag != 10)
  {
    cout << "Non-fatal problem with FBA! " << exitflag << std::endl;
  }
}

void Cell::do_fba(QVector<BactParam> bps)
{
  /* This method performs flux balance analysis. It solves the linear programming problem:
   *  S*f=x, where S is an m by n stoichiometric matrix, f the flux vector (n by 1)(also exchange fluxes!) and x=dm/dt (the change of metabolite concentration over time, an m by 1 vecor)
   *  with the constraints: lb(i)<f(i)<ub(i) and x(i)=0 for all internal and external metabolites, and unconstrained for exchange metabolites
   *  (to account for possible influx and efflux in or out of the system). The function obj*f (where f is a n by 1 objective vector) is optimized given the constraints. The stoichiometrix matrix
   *  is given in sparse form, using three vectors iRow, jCol and sparseS (see glpk manual).
   */

  int exitflagtemp;

  // define problem

  lp = bps.at(indexnr - 1).lp;
  params = bps.at(indexnr - 1).params;

  glp_set_prob_name(lp, "linprog");
  glp_set_obj_dir(lp, GLP_MAX); // maximize objective function
  // glp_term_out(GLP_OFF);//suppress output of glpk

  // Set reaction bounds: this changes during simulation, so has to be done here.
  for (int i = 1; i <= numreactions; i++)
  {
    double ubtemp = ub[i - 1];
    double lbtemp = lb[i - 1];
    if (ub[i - 1] < lb[i - 1])
      std::cout << "Problem in GLPK: ub<lb: " << lbtemp << " " << ubtemp << endl;
    else
    {
      ubtemp = min(100 * par.max_metup, ubtemp);
      lbtemp = min(100 * par.max_metup, lbtemp);
      if (ubtemp < 0.000001)
        ubtemp = 0;
      if (lbtemp < 0.000001)
        lbtemp = 0;
      if (ubtemp > lbtemp)
        glp_set_col_bnds(lp, i, GLP_DB, lbtemp, ubtemp);
      else
        glp_set_col_bnds(lp, i, GLP_FX, lbtemp, ubtemp);
    }

    if (obj[i - 1] != 0)
    {
      glp_set_obj_coef(lp, i, obj[i - 1]);
    }
  }

  // glp_set_row_bnds(lp, nummetabolites+1, GLP_DB, 0, maxflux);

  exitflagtemp = glp_simplex(lp, params);         // perform simplex method
  growthrate = glp_get_col_prim(lp, bm_reaction); // assign growth rate
  double originalgrowthrate = growthrate;
  growthrate = growthrate / par.growthmod;

  for (int i = 0; i < numreactions; i++)
  {
    flux[i] = glp_get_col_prim(lp, i + 1); // assign flux values
    if (used[i] == 0 && flux[i] > 0.0000000000000001)
      used[i] = 1;
  }

  double carbonflux = 0;
  double hydrogenflux = 0;
  double nitrogenflux = 0;
  double oxygenflux = 0;
  energyflux = 0; // This is a class variable
  double energyinflux = 0;
  double energyoutflux = 0;
  Exchange exch;
  for (auto const &x : exchange_dict)
  {
    exch = x.second;
    if (exch.reaction_in != 0)
    {
      double influx = flux[exch.reaction_in - 1];
      double outflux = flux[exch.reaction_out - 1];
      carbonflux -= influx * (exch.carbon);
      carbonflux += outflux * (exch.carbon);
      hydrogenflux -= influx * (exch.hydrogen);
      hydrogenflux += outflux * (exch.hydrogen);
      nitrogenflux -= influx * (exch.nitrogen);
      nitrogenflux += outflux * (exch.nitrogen);
      oxygenflux -= influx * (exch.oxygen);
      oxygenflux += outflux * (exch.oxygen);
      energyinflux += influx * (exch.energy);
      energyoutflux += outflux * (exch.energy);
      if (exch.energy == 0.0)
      {
        // if(x.second != "h"){
        if (influx > 0.00001 || outflux > 0.00001)
        {
          if (x.first != "h")
          {
            std::cout << x.first << " not found in energy file, used here at " << influx << " flux" << std::endl; // Edited out when using mucus, because we don't have energy values for that
          }
        }
        //}
      }
    }
  }

  // Compensate for mass lost in creating biomass
  nitrogenflux += originalgrowthrate * exchange_dict[bps.at(indexnr - 1).speciesName].nitrogen;
  carbonflux += originalgrowthrate * exchange_dict[bps.at(indexnr - 1).speciesName].carbon;
  hydrogenflux += originalgrowthrate * exchange_dict[bps.at(indexnr - 1).speciesName].hydrogen;
  oxygenflux += originalgrowthrate * exchange_dict[bps.at(indexnr - 1).speciesName].oxygen;

  bool displayexchanges = false; // by default. If this turns true, we will display the exchanges
  if (par.alwayssaveenergy)
  {
    displayexchanges = true;
  }
  energyflux = energyoutflux - energyinflux; // Should be negative
  if (abs(carbonflux) > 0.00001)
  {
    std::cout << "net carbon flux is " << carbonflux << " umol for " << species << std::endl;
    std::cout << "growth of " << originalgrowthrate << " for " << species << std::endl;
    displayexchanges = true;
  }
  if (abs(hydrogenflux) > 0.00001)
  {
    std::cout << "net hydrogen flux is " << hydrogenflux << " umol for " << species << std::endl;
    std::cout << "growth of " << originalgrowthrate << " for " << species << std::endl;
    displayexchanges = true;
  }
  if (abs(nitrogenflux) > 0.00001)
  {
    std::cout << "net nitrogen flux is " << nitrogenflux << " umol for " << species << std::endl;
    std::cout << "growth of " << originalgrowthrate << " for " << species << std::endl;
    displayexchanges = true;
  }
  if (abs(oxygenflux) > 0.00001)
  {
    std::cout << "net oxygen flux is " << oxygenflux << " umol for " << species << std::endl;
    std::cout << "growth of " << originalgrowthrate << " for " << species << std::endl;
    displayexchanges = true;
  }
  if ((energyoutflux > energyinflux) && growthrate > 0.00001)
  {
    std::cout << "net energy flux is " << energyflux << " mjoule for " << species << std::endl;
    std::cout << "growth of " << growthrate << " for " << species << std::endl;
    displayexchanges = true;
  }

  if (displayexchanges)
  {
    for (auto const &x : exchange_dict)
    {
      exch = x.second;
      if (exch.reaction_in != 0)
      {
        if ((flux[exch.reaction_in - 1] > 0.00000001) && (flux[exch.reaction_in - 1] > flux[exch.reaction_out - 1]))
        {
          std::cout << x.first << " " << flux[exch.reaction_in - 1] - flux[exch.reaction_out - 1] << " in at " << exch.energy << " energy" << std::endl;
        }
        if ((flux[exch.reaction_out - 1] > 0.00000001) && (flux[exch.reaction_out - 1] > flux[exch.reaction_in - 1]))
        {
          std::cout << x.first << " " << flux[exch.reaction_out - 1] - flux[exch.reaction_in - 1] << " out at " << exch.energy << " energy" << std::endl;
        }
      }
    }
  }

  // write energies to file (SLOW!)
  bool writeenergies = displayexchanges; // only write energies when there's a problem
  // displayexchanges = false;
  if (writeenergies)
  {
    if (growthrate > par.MIN_GROWTH)
    {
      ofstream energyfile;
      energyfile.open(simulationdir + "/energyperpopperstep.txt", fstream::app);
      energyfile << energyinflux << "," << energyoutflux << "," << originalgrowthrate << std::endl;
      energyfile.close();
    }
    if (displayexchanges)
    {
      ofstream imbalancefile;
      imbalancefile.open(simulationdir + "/imbalancedreactions.txt", fstream::app);
      imbalancefile << "Imbalance of " << energyflux << " for " << species << " at " << growthrate << " growth, " << energyinflux << " influx" << energyoutflux << " outflux" << std::endl;
      imbalancefile << "Imbalance of " << carbonflux << " carbon " << hydrogenflux << " hydrogen " << nitrogenflux << " nitrogen " << oxygenflux << " oxygen" << std::endl;
      for (auto const &x : react_dict)
      {
        if (used[x.second - 1] == 1)
        {
          if (flux[x.second - 1] > 0.0000000000000001)
          {
            imbalancefile << x.first << "," << flux[x.second - 1] << std::endl;
          }
        }
      }
      imbalancefile << std::endl;
      imbalancefile.close();
    }
  }

  // Can be used to print the exact reactions to console when flux is unbalanced
  /*
    for (auto const& x : exchange_dict)
    {
        exch = x.second;
        if (exch.reaction_in != 0)
        {
        if(flux[exch.reaction_in - 1]>0.0001){
        std::cout<<x.first<<" "<<flux[exch.reaction_in - 1]<<" at "<< exch.nitrogen<<" carbon equals " << flux[exch.reaction_in - 1] * (double)(exch.nitrogen)<<" in"<<std::endl;
        }
        if(flux[exch.reaction_out - 1]>0.0001){
        std::cout<<x.first<<" "<<flux[exch.reaction_out - 1]<<" equals "<< exch.nitrogen<<" carbon equals "<<flux[exch.reaction_out - 1] * (double)(exch.nitrogen)<<" carbon out"<<std::endl;
        }
        }
    }

    std::cout <<"And biomass contributes " << originalgrowthrate <<" times " <<exchange_dict[bps.at(indexnr-1).speciesName].nitrogen <<" equals "<<originalgrowthrate * exchange_dict[bps.at(indexnr-1).speciesName].nitrogen<<std::endl;
    */
  if (growthrate < 0)
  {
    ofstream errorfile;
    errorfile.open(simulationdir + "/errorlog.txt", fstream::app);
    errorfile << "Negativegrowtherror"
              << "," << species << "," << growthrate << std::endl;
    errorfile.close();
    growthrate = 0;
  }
  switch (exitflagtemp) // Determine if an error has occurred during the FBA.
  {
  case 0:
  {
    if (displayexchanges)
    {
      exitflag = -1; // succesful solution, but not thermodynamically correct and/or has a mass imbalance
    }
    exitflag = 0; // The LP problem instance has been successfully solved, no output
    break;
  }
  case GLP_EBADB:

  {
    std::cout << "Unable to start the search, because the initial basis specified in the problem object is invalid.\n";
    exitflag = 1;
    break;
  }
  case GLP_ESING:
  {
    std::cout << "Unable to start the search for " << species << ", because the basis matrix corresponding to the initial basis is singular within the working precision.\n";
    std::cout << "This run will now abort";
    exitflag = 2;
    break;
  }
  case GLP_ECOND:
  {
    std::cout << "Unable to start the search, because the basis matrix corresponding to the initial basis is ill-conditioned, i.e. its condition number is too large.\n";
    exitflag = 3;
    break;
  }
  case GLP_EBOUND:
  {
    std::cout << "Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds.\n";
    exitflag = 4;
    break;
  }
  case GLP_EFAIL:
  {
    std::cout << "The search was prematurely terminated due to the solver failure.\n";
    exitflag = 5;
    break;
  }
  case GLP_EOBJLL:
  {
    std::cout << "The search was prematurely terminated, because the objective function being maximized has reached its lower limit and continues decreasing.\n";
    exitflag = 6;
    break;
  }
  case GLP_EOBJUL:
  {
    std::cout << "The search was prematurely terminated, because the objective function being minimized has reached its upper limit and continues increasing.\n";
    exitflag = 7;
    break;
  }
  case GLP_EITLIM:
  {
    std::cout << "The search was prematurely terminated, because the simplex iteration limit has been exceeded.\n";
    exitflag = 8;
    break;
  }
  case GLP_ETMLIM:
  {
    std::cout << "The search for " << species << " was prematurely terminated, because the time limit has been exceeded.\n";
    exitflag = 9;
    break;
  }
  case GLP_ENOPFS:
  {
    std::cout << "The LP problem instance has no primal feasible solution.\n";
    exitflag = 10;
    break;
  }
  case GLP_ENODFS:
  {
    std::cout << "The LP problem instance has no dual feasible solution\n.";
    exitflag = 11;
    break;
  }
  }

  assert(exitflag != 2); // Crash if it is, run will be ruined by slowdown anyway if we let it continue
}

void Cell::assignUB(std::string key, double value)
{
  /* This method can be used to set the upper bound of an in reaction for a metabolite (key) to a specific (value). */
  std::map<std::string, Exchange>::iterator it = exchange_dict.find(key);
  Exchange exch = it->second;
  if (exch.reaction_in != 0 && it != exchange_dict.end())
    ub[exch.reaction_in - 1] = value;
}

void Cell::assignInFlux(double &var, std::string key)
{
  /* This method can be used to assign the flux of an in reaction for a metabolite (key) to a variable (&var). */
  std::map<std::string, Exchange>::iterator it = exchange_dict.find(key);
  Exchange exch = it->second;
  if (exch.reaction_in != 0 && it != exchange_dict.end())
    var = flux[exch.reaction_in - 1];
}

void Cell::assignOutFlux(double &var, std::string key)
{
  /* This method can be used to assign the flux of an out reaction for a metabolite (key) to a variable (&var). */
  std::map<std::string, Exchange>::iterator it = exchange_dict.find(key);
  Exchange exch = it->second;
  if (exch.reaction_out != 0 && it != exchange_dict.end())
    var = flux[exch.reaction_out - 1];
}

void Cell::Init(int i2, int j2)
{
  /* This method initializes the Cell, no longer reads in the metabolic network as a matrix to be used by GLPK.
   * It also contains a lot of test on how fast the bacterium can grow on different metabolites that can be used for debugging and such.
   */

  // Metabolic network initialized

  for (int k = 0; k < numreactions; k++)
  {
    flux[k] = 0;
    fluxtot[k] = 0;
  }

  // std::cout<<"Initialized!,";

  ipos = i2; // sets x-pos of cell
  jpos = j2; // sets y-pos of cell
  numancestors = 0;
  growthrate = 0;
  volume = par.initialvolume;
  exist = true;
  hasmoved = false;
  active = true;
  exitflag = 0;

  ancestorlabel = ancestorcount;
  firstancestorlabel = ancestorlabel;
}

void Cell::Die()
{
  /* This method lets the bacterium die and makes the position empty. */
  Cell Celltemp;
  *this = Celltemp;
}

double Cell::Crossfeeding()
{
  /* This method calculates how much carbon bacteria have aquired through cross feeding. */
  double crossflux = 0;
  Exchange exch;
  for (auto const &x : exchange_dict)
  {
    if (x.first.compare("glc_D") != 0 && x.first.compare("fru") != 0)
    {
      exch = x.second;
      if (exch.reaction_in != 0 && exch.reaction_out != 0)
      {
        crossflux += exch.carbon * max(fluxtot[exch.reaction_in - 1] - fluxtot[exch.reaction_out - 1], 0) / (double)(par.timescale);
      }
      else if (exch.reaction_in != 0)
      {
        crossflux += exch.carbon * fluxtot[exch.reaction_in - 1] / (double)(par.timescale);
      }
    }
  }
  return crossflux;
}

double Cell::Totflux()
{
  /* This method calculates the total amount of carbon the bacteria have aquired. */
  double totflux = 0;
  Exchange exch;
  for (auto const &x : exchange_dict)
  {
    exch = x.second;
    if (exch.reaction_in != 0)
    {
      totflux += exch.carbon * fluxtot[exch.reaction_in - 1] / (double)(par.timescale);
    }
  }
  return totflux;
}

void Cell::ShowInfo()
{
  /* This method can be used to print some info on a bacterium such as the growth rates on different metabolites, used for debugging. */
  std::cout << "The species is: " << species << std::endl;
  // for (int i = 0; i < ub.size(); ++i) std::cout << ub.at(i) << " ";
  // std::cout << std::endl;
  // for (auto const& x : exchange_dict) std::cout << x.first << " in reaction: " << x.second.reaction_in << " out reaction: " << x.second.reaction_out << std::endl;
  // std::cout << std::endl;
  // for (int i = 0; i < flux.size(); ++i) std::cout << flux.at(i) << " ";
  // std::cout << std::endl;
  for (auto const &x : react_dict)
    if (used[x.second - 1] == 1)
      std::cout << x.first << ", flux: " << flux[x.second - 1] << std::endl;
  // for (int i = 0; i < nummetabolites; ++i) std::cout << glp_get_row_prim(lp, i+1) << " ";
  std::cout << std::endl;
}
