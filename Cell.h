#ifndef CELL_H
#define CELL_H

#include <QVector>
#include "project.h"
#include <string>
#include "BactParam.h"

/* This class implements the bacterium (metapopulation of bacteria) that occupies one position.
 * The class is constructed on the basis of an instance of class BactParam that contains the metabolic network for a specific species.
 * This class contains methods to use that metabolic network through FBA to grow and consume and produce metabolites.
 */

class Cell
{
    friend class Grid;
    friend class Gut_Output;
    friend class Graphics;

private:
    QVector<int> genome;     // vector of length numreactions with 1 if a reaction is on (now just 1 for all, so not really used)
    QVector<int> used;       // vector of length numreactions with 1 if a reaction is used and 0 if not (mostly for debugging)
    QVector<double> obj;     // vector of length numreactions with 1 for the biomass reaction and 0 for all others (used in FBA)
    QVector<double> SparseS; // Stoichiometric values (forms the sparse matrix S together with iRow and jCol)
    QVector<int> iRow;       // Row index (metabolite) (forms the sparse matrix S together with SparseS and jCol)
    QVector<int> jCol;       // Column index (reaction) (forms the sparse matrix S together with iRow and SparseS)
    QVector<double> lb;      // vector of length numreactions with lower bounds of reactions (used in FBA)
    QVector<double> ub;      // vector of length numreactions with upper bounds of reactions (used in FBA)
    QVector<double> lbubfood;
    QVector<double> flux;    // vector of length numreactions with fluxes of reactions (used in FBA)
    QVector<double> fluxtot; // vector of length numreactions with summation of fluxes over lifetime of bacterium

    std::map<std::string, Exchange> exchange_dict; // map: key is the id of the metabolite, value is number of reactions in, number of reactions out and number of carbon atoms in molecule.
    QVector<int> exchange_list;                    // vector with exchange metabolite ids

    std::map<std::string, int> react_dict; // map: key is reaction id, value is index of reaction

    glp_prob *lp;
    glp_smcp *params;

    std::string species; // species name, taken from SBML model
    int indexnr;         // corresponding to the order in which the bacteria are added (index in b_params in simulation.cpp)
    int numreactions;    // number of reactions
    int nummetabolites;  // number of metabolites
    int numstoich;       // number of stoichiometric values (non-zero values in matrix S)
    double maxflux;      // not used
    int bm_reaction;     // the index of the biomass reaction in vector flux and similar vectors

    double growthrate; // the rate by which the cell grows, the flux of the biomass reaction
    double energyflux; // The calculated thermodynamic balance associated with the current fluxes

    // some values calculated in the init used for debugging
    double volume; // The volume of the bacterium/metapopulation
    bool exist;    // whether there is a living bacterium in Cell or just an empty position
    bool hasmoved; // whether the bacterium has just moved
    bool active;
    int exitflag;           // exit code of GLPK
    int ancestorlabel;      // not used
    int firstancestorlabel; // not used
    int ipos;               // x position in the grid
    int jpos;               // y position in the grid
    int timebirth;          // iteration at which cell is created, not used
    int time_alive;         // iterations that cell has been alived, not used
    int ancestor;           // direct ancestor
    int genomesize;         // number of reactions of bacterium, not used
    int genomesizeinit;     // number of reactions of bacterium at intialization, not used

public:
    int numancestors;               // number of ancestors
    Cell();                         // Method to construct an empty cell without a bacterium, so an empty position.
    Cell(BactParam bp);             // Method to construct a Cell with a bacterium, the species is determined by the argument bp, which contains the metabolic network.
    Cell(const Cell &c);            // Method to construct a Cell on the basis of another instance of Cell.
    ~Cell();                        // Destructor of Cell class.
    Cell &operator=(const Cell &c); // Method to create a reference to the Cell class.
    double oxygenation;
    bool CanGrow() const;                             // Determines whether FBA is performed, to speed-up simulations, not used!!
    void Do_Metabolism(QVector<BactParam> bps);       // This method calls the method for FBA, sees if a threshold for minimal growth is achieved and calculates total fluxes.
    void do_fba(QVector<BactParam> bps);              // This method performs flux balance analysis. It solves the linear programming problem: S*f=x.
    void assignUB(std::string key, double value);     // This method can be used to set the upper bound of an in reaction for a metabolite (key) to a specific (value).
    void assignInFlux(double &var, std::string key);  // This method can be used to assign the flux of an in reaction for a metabolite (key) to a variable (&var).
    void assignOutFlux(double &var, std::string key); // This method can be used to assign the flux of an out reaction for a metabolite (key) to a variable (&var).
    void Init(const int i, const int j);              // This method initializes the Cell, specifically reads in the metabolic network as a matrix to be used by GLPK. It also contains a lot of test on how fast the bacterium can grow on different metabolites that can be used for debugging and such.
    void Die();                                       // This method lets the bacterium die and makes the position empty.
    double Crossfeeding();                            // This method calculates how much carbon bacteria have aquired through cross feeding.
    double Totflux();                                 // This method calculates the total amount of carbon the bacteria have aquired.
    void ShowInfo();                                  // This method can be used to print some info on a bacterium such as the growth rates on different metabolites, used for debugging.
};

#endif
