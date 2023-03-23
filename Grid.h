#ifndef GRID_H
#define GRID_H

#include <map>
#include <set>
#include <string>
#include <QVector>
#include <fstream>
#include "Cell.h"
#include "BactParam.h"

/** This class implements the grid which contains the bacteria and metabolite concentrations.
 */

class Grid
{
    friend class Graphics;
    friend class Gut_Output;

private:
    int nrow;                                       // Number of rows for the grid
    int ncol;                                       // Number of columns for the grid
    int numexchangemet;                             // Number of metabolites in the system (all metabolites present in the SBML models).
    std::map<std::string, int> exchange_metab_dict; // the names of all the metabolites in the system mapped to a index (indexes are easier in for loops and such)

    Cell ***Cellfield;  // contains the bacteria
    double ***metfield; // contains the concentrations of metabolites
    std::map<std::string, double> ***enzymefield;
    double *sumconc; // This array contains the average concentrations of the metabolites through the whole system.
    void Record_Host_Uptake(int metabolite, double uptake);
    void Record_Fecal_Diffusion(int metabolite, double uptake);

public:
    std::map<std::string, double> fecesbacmap;                 // Record all bacteria excreted per step
    std::map<std::string, double> fecesmetmap;                 // Record all metabolites excreted per step
    std::map<std::string, double> inputmetmap;                 // Record all metabolites entered into the system
    std::map<int, double> uptakemetmap;                        // Record all metabolites taken up per step
    std::vector<std::string> bacspecies;                       // Hold all species used in bacplace
    Grid(int nr, int nc, QVector<std::string> exchange_metab); // This method constructs the grid with dimension nr X nc.
    Grid(const Grid &g);                                       // Method that creates a Grid class on the basis of another instance of the Grid class.
    ~Grid();                                                   // Destructor of the Grid class.
    Grid &operator=(const Grid &c);

    void Init(QVector<BactParam> bps); // Initializer for the Grid class. This method fills the Cellfield with bacteria according to a certain initial density and with an equal distribution of the different species.

    double CountCellNeighbour(int i, int j) const;               // This method counts how many living neighbours a bacteria at a certain position has.
    double *MetFieldNeighbour(int i, int j, int k, int n) const; // Get concentration of a certain metabolite (k) on a position neighbouring position (i, j) in a certain direction (n).
    void Update_Sumconc();                                       // This method calculates the concentrations of the different metabolites averaged over the entire grid.
    double Avggrowth() const;                                    // This method calculates the average growth rate of all bacteria.
    int Numbercells() const;                                     // This method calculates the number of bacteria alive in the grid.

    void Diffuse_Cells();                       // This method regulates the movement of the bacteria.
    void Diffuse_Metabolites(double diffconst); // This method calculates the diffusion of metabolites through the grid. At the borders it also determines the uptake by the host.
    void DriftCells(QVector<BactParam> bps);    // This method regulates the drift of the metabolites flowing through the system and if applicable also the drift of the bacteria.
    void CellShuffle();                         // This method randomly shuffles the bacteria with a position anywhere in the grid. Normally not used!
    void Mix_Metabolites();                     // This method sets the concentrations of all metabolites equal to the average. Normally not used!
    void ReadEnv(int i, int j);                 // This method sets the upper bounds of the uptake reactions of the bacteria according to the locally available concentrations of metabolites.
    void Update_Concentrations(int i, int j);   // This method updates the local concentrations according to what the bacteria at that position take up and excrete.
    void Place_Bacterium(QVector<BactParam> bps);

    void Set_Unlimited_Metabolite();                                         // This method sets the concentrations of certain metabolites, namely H+, H2, h20, na1, nh4, pi and so4, to high values.
    void MixedSugarInflux();                                                 // This method adds sugar to the left end of the grid, how much is determined by the parameter foodinglucose and the ratio between glucose and fructose by frac_fru.
    void HostInflux();                                                       // This method adds oxygen at the positions at the borders. How much oxygen is determined by the parameter o2.
    bool Cellexist(int i, int j, int k) { return Cellfield[i][j][k].exist; } // Determine if position is occupied by living bacterium.
    void Step_Metabolism(QVector<BactParam> bps);                            // Let all bacteria in the entire grid metabolize and update the concentrations of metabolites accordingly.
    void Grid_Die_Reproduce();                                               // This method regulates density dependent cell death and division of the bacteria.
    void UsedReactions();                                                    // Print which reactions of their metabolic network the bacteria in the grid actually use. Used for debugging and such.
};
#endif
