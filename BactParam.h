#ifndef BACTPARAM_H
#define BACTPARAM_H

#include <QVector>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <pngwriter.h>
#include <sbml/SBMLTypes.h>
#include "project.h"

// This struct can contain the index of the in reaction, the out reaction and the number of carbon atoms for a metabolite.
struct Exchange
{
    int reaction_in;
    int reaction_out;
    double carbon;
    double hydrogen;
    double nitrogen;
    double oxygen;
    double energy;
};

// functions implemented in BP_functions.cpp
bool hasPrefix(std::string const &fullString, std::string const &prefix);
bool hasEnding(std::string const &fullString, std::string const &ending);
std::string make_universal_id(std::string id);
int getCarbon(std::string notes);
int getHydrogen(std::string notes);
int getNitrogen(std::string notes);
int getOxygen(std::string notes);
std::string getFullName(std::string notes);
std::string getKeggName(std::string notes);
std::string getChebiName(std::string notes);
std::string getName(Model *m, std::string id);
int contains_pair(QVector<int> v1, QVector<int> v2, int x1, int x2);

/* The class BactParam contains methods to read in the SBML model.
 * The information of the SBML model is saved in a universal way such that the information can be used in the class Cell.
 * This class is passed as a parameter to class Cell.
 */
class BactParam
{
public:
    int indexnr;             // corresponding to the order in which the bacteria are added (index in b_params in simulation.cpp)
    double FW;               // value for uniform flux weights (same value for all reactions)
    std::string speciesName; // species name, taken from SBML model
    int bm_reaction;         // the index of the biomass reaction in vector flux and similar vectors
    QVector<int> genome;     // vector of length numreactions with 1 if a reaction is on (now just 1 for all, so not really used)
    QVector<int> used;       // vector of length numreactions with 1 if a reaction is used and 0 if not (mostly for debugging)
    QVector<double> lbs;     // vector of length numreactions with lower bounds of reactions (used in FBA)
    QVector<double> ubs;     // vector of length numreactions with upper bounds of reactions (used in FBA)
    QVector<double> fluxs;   // vector of length numreactions with fluxes of reactions (used in FBA)
    QVector<double> objs;    // vector of length numreactions with 1 for the biomass reaction and 0 for all others (used in FBA)
    QVector<int> iRow;       // Row index (metabolite) (forms the sparse matrix S together with SparseS and jCol)
    QVector<int> jCol;       // Column index (reaction) (forms the sparse matrix S together with iRow and SparseS)
    QVector<double> SparseS; // Stoichiometric values (forms the sparse matrix S together with iRow and jCol)z
    glp_prob *lp;
    glp_smcp *params;
    std::map<std::string, int> metab_dict;         // map: key is metabolite id, value is index of metabolite
    std::map<std::string, int> react_dict;         // map: key is reaction id, value is index of reaction
    int numreactions;                              // number of reactions
    std::map<std::string, Exchange> exchange_dict; // map: key is the id of the metabolite, value is an instance of struct Exchange.
    QVector<int> exchange_list;                    // vector with indexes of exchange reactions
    int nummetabolites;                            // number of metabolites
    int numstoich;                                 // number of stoichiometric values (non-zero values in matrix S)
};

/* The classes that inherit from BactParam are largely similar, but their methods are adapted to a specific SBML model.
 * The ids of the reactions and metabolites of the models they read in will be made universal, so that the ids of the same metabolites in different models correspond to each other.
 */

class AGORAParam : public BactParam
{
public:
    AGORAParam(std::string file, int i, double fw, std::map<std::string, double> gibbsmap); // Constructs a EcoliParam object on the basis of an SBML file (file) and a value for the flux weight (fw).
    ~AGORAParam();                                                                          // Destructor, does nothing in this case.
    void remove_reaction(std::string);                                                      // Can be used to remove any desired reaction

private:
    void createDicts(Model *m, std::map<std::string, double> gibbsmap); // This method fills the maps/dictionairys react_dict, metab_dict and exchange_dict.
    void addReaction(Reaction *react);                                  // This method adds the reactions of react_dict to the vectors genome, lbs, ubs, fluxs and objs.
    std::vector<double> addReactants(Model *m, Reaction *react);        // This method adds the stoichiometric values of the reactants (<0) to the matrix S.
    std::vector<double> addProducts(Model *m, Reaction *react);         // This method adds the stoichiometric values of the products (>0) to the matrix S.
    void setBiomassReaction(Model *m);                                  // This method creates the standard biomass reaction used to replace the biomass reaction present in the SBML model.

    void add_FWs();              // This method adds the flux weights (one value for all reactions in this case) to the matrix S.
    void addBacterium(Model *m); // This method calls all the methods to load one bacterium species into the class.
    std::ofstream EXCHANGENAMEFILE;
    std::ofstream EXCHANGEIDFILE;
    std::ofstream EXCHANGEKEGGFILE;
    std::ofstream EXCHANGECHEBIFILE;
};

#endif
