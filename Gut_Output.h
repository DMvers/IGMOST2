#ifndef FILE_OUTPUT_H
#define FILE_OUTPUT_H

#include <fstream>
#include "Grid.h"
#include <sys/types.h>
#include <sys/stat.h>

/* This class takes care of the output files in the data folder. */

class Gut_Output
{

private:
    std::ofstream EXCHANGEFILE;  // handle of the output file "exchange.dat"
    std::ofstream SCFAFILE;      // handle of the output file "scfa.csv"
    std::ofstream BACTERIAFILE;  // handle of the output file "bacteriasummary"
    std::ofstream METABFILE;     // handle of the output file "metabolites"
    std::ofstream PARAMFILE;     // handle of the output file "parameters"
    std::ofstream REACTFILE;     // handle of the output file "reactions"
    std::ofstream FECESBACFILE;  // handle of the output file "fecesbac"
    std::ofstream FECESMETFILE;  // handle of the output file "fecesmet"
    std::ofstream UPTAKEMETFILE; // handle of the output file "uptake"
    std::ofstream TOTENERGYFILE;
    std::ofstream TOTELEMENTFILE;

    std::string datadir; // Name of the data folder

public:
    std::map<std::string, double> bacteriamap;
    ~Gut_Output();                       // Destructor of the Gut_Output class closes the handles of the files.
    void Init();                         // This method initializes the class Gut_Output.
    void Write_Output(Grid &g);          // This method writes to the output file "exchange.dat"
    void Write_SCFA(Grid &g);            // This method writes to the output file "scfa.csv"
    void Write_BacteriaSummary(Grid &g); // Print bacteria count per species (also present in exchange.dat, but easier to read in this format)
    void Write_Metabolites(Grid &g);     // Print which metabolites are present
    void Write_Params();                 // Print some important parameters, including the random seed
    void Write_Reactions(Grid &g);       // Print the reactions used by each population in the last step
    void Write_Feces(Grid &g);           // Print outflow over the last period
    void Write_Host_Uptake(Grid &g);     // Print host uptake (normally not present)
    void Write_Energy(Grid &g, std::map<std::string, double> gibbsmap);
    void Write_Elements(Grid &g, std::map<std::string, double> carbonmap, std::map<std::string, double> hydrogenmap, std::map<std::string, double> nitrogenmap, std::map<std::string, double> oxygenmap);
};
#endif
