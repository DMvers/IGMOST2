#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include "Gut_Output.h"
#include "parameter.h"

extern Par par;
extern char *SpeciesName;
extern int simultime;
extern int numexchangemet;
extern int exchangeint;
extern std::string simulationdir;

/* This class takes care of the output files in the data folder. */

void Gut_Output::Init()
{
    /* This method initializes the class Gut_Output.
     * In this method the directory data is created with
     * the empty files "exchange.dat" and "scfa.csv" in it.
     */
    datadir = simulationdir;
    datadir += par.datadir;
    // strcpy(datadir,simulationdir);
    // strcat(datadir,par.datadir);

    // workaround because mkdir is actually deprecated
    char ddir[datadir.size() + 1];
    strcpy(ddir, datadir.c_str());
    if (mkdir(ddir, 0777) == -1) // creating a directory
    {
        cerr << "Directory " << datadir << " already exists" << endl;
    }

    std::string exchangefilename;
    std::string scfafilename;
    std::string bacteriafilename;
    std::string metabfilename;
    std::string paramfilename;
    std::string reactionfilename;
    std::string fecesbacfilename;
    std::string fecesmetfilename;
    std::string uptakemetfilename;
    std::string energyfilename;
    std::string elementfilename;
    // std::cout << "data directory = " <<datadir <<std::endl;
    // std::cout << "sim directory = " <<simulationdir <<std::endl;

    // Remove prefixes
    int startofname = par.testname.find_last_of('/');
    std::string onlyname = par.testname.substr(startofname + 1, par.testname.size()); // This still works if there is no '/', returns entire string

    exchangefilename = datadir + "/exchange" + onlyname + to_string(par.runid) + ".dat";
    scfafilename = datadir + "/scfa" + onlyname + to_string(par.runid) + ".csv";
    bacteriafilename = simulationdir + "/bacteriasummary" + onlyname + to_string(par.runid) + ".csv";
    metabfilename = datadir + "/metabolites" + onlyname + to_string(par.runid) + ".txt";
    paramfilename = datadir + "/parameters" + onlyname + to_string(par.runid) + ".txt";
    reactionfilename = datadir + "/reactions" + onlyname + to_string(par.runid) + ".txt";
    fecesbacfilename = datadir + "/fecesbac" + onlyname + to_string(par.runid) + ".txt";
    fecesmetfilename = datadir + "/fecesmet" + onlyname + to_string(par.runid) + ".txt";
    uptakemetfilename = datadir + "/uptake" + onlyname + to_string(par.runid) + ".txt";
    energyfilename = datadir + "/energy" + onlyname + to_string(par.runid) + ".txt";
    elementfilename = datadir + "/element" + onlyname + to_string(par.runid) + ".txt";

    std::cout << "param directory = " << paramfilename << std::endl;

    // Open files
    EXCHANGEFILE.open(exchangefilename);
    SCFAFILE.open(scfafilename);
    BACTERIAFILE.open(bacteriafilename);
    METABFILE.open(metabfilename);
    PARAMFILE.open(paramfilename);
    REACTFILE.open(reactionfilename);
    FECESMETFILE.open(fecesmetfilename);
    FECESBACFILE.open(fecesbacfilename);
    UPTAKEMETFILE.open(uptakemetfilename);
    TOTENERGYFILE.open(energyfilename);
    TOTELEMENTFILE.open(elementfilename);

    // Write some headers
    EXCHANGEFILE << "timestep;x;y;species;volume;growth;energy";
    for (std::string m : par.outputmetabs)
    {
        EXCHANGEFILE << ";" << m;
    }
    EXCHANGEFILE << endl;
    TOTENERGYFILE << "timestep,in,system,out" << std::endl;
    TOTELEMENTFILE << "timestep,cin,csystem,cout,hin,hsystem,hout,nin,nsystem,nout,oin,osystem,oout" << std::endl;
}

Gut_Output::~Gut_Output()
{
    /* Destructor of the Gut_Output class closes the handles of the files. */
    EXCHANGEFILE.close();
    SCFAFILE.close();
    BACTERIAFILE.close();
    METABFILE.close();
    PARAMFILE.close();
    REACTFILE.close();
    FECESMETFILE.close();
    FECESBACFILE.close();
    UPTAKEMETFILE.close();
    TOTENERGYFILE.close();
    TOTELEMENTFILE.close();
}

void Gut_Output::Write_Output(Grid &g)
{
    /* This method writes to the output file "exchange.dat"
     * with data about the bacteria and their exchange fluxes.
     */
    for (int i = 1; i <= g.nrow; i++)
    {
        for (int j = 1; j <= g.ncol; j++)
        {
            for (int k2 = 0; k2 < par.maxcells; k2++)
            {
                if (g.Cellfield[i][j][k2].exist == true)
                {
                    // double cross;
                    // cross = (g.Cellfield[i][j][k2].Totflux() > 0.1) ? g.Cellfield[i][j][k2].Crossfeeding() / g.Cellfield[i][j][k2].Totflux() : 0;
                    EXCHANGEFILE << simultime << ";" << i << ";" << j << ";" << g.Cellfield[i][j][k2].species << ";" << g.Cellfield[i][j][k2].volume << ";" << g.Cellfield[i][j][k2].growthrate << ";" << g.Cellfield[i][j][k2].energyflux;
                    for (std::string m : par.outputmetabs)
                    {
                        if (g.Cellfield[i][j][k2].exchange_dict.count(m) > 0)
                        {
                            int r_in = g.Cellfield[i][j][k2].exchange_dict.at(m).reaction_in;
                            int r_out = g.Cellfield[i][j][k2].exchange_dict.at(m).reaction_out;
                            double s = g.Cellfield[i][j][k2].volume * par.cellfrac / par.timescale;
                            double influx = (r_in != 0) ? g.Cellfield[i][j][k2].flux.at(r_in - 1) * s : 0;
                            double outflux = (r_out != 0) ? g.Cellfield[i][j][k2].flux.at(r_out - 1) * s : 0;
                            if (influx > 0 && outflux > 0)
                            {
                                // std:cout << "something's weird with the flux " <<influx<<outflux<<m;
                            }
                            double netflux = influx - outflux;
                            if (((netflux < par.DIFF_CUTOFF) && (netflux > 0)) || ((netflux > -1 * par.DIFF_CUTOFF) && netflux < 0))
                            {
                                netflux = 0;
                            }
                            EXCHANGEFILE << ";" << (netflux);
                            // EXCHANGEFILE <<";"<< m << ";" << influx << ";" << outflux;
                        }
                        else
                        {
                            EXCHANGEFILE << ";0";
                        }
                    }
                    EXCHANGEFILE << endl;
                }
            }
        }
    }
}

void Gut_Output::Write_SCFA(Grid &g)
{
    /* This method writes to the output file "scfa.csv"
     * with the concentrations of the SCFAs at every position.
     */
    for (int j = 1; j <= g.ncol; j++)
    {
        for (int i = 0; i <= g.nrow + 1; i++)
        {
            SCFAFILE << simultime << "," << i << "," << j;
            for (std::string m : {"lcts", "2FuLa", "lac_L", "lac_D", "ac","12ppd_S","but"}) //,"for", "succ"})
            {
                double metab = g.metfield[i][j][g.exchange_metab_dict.at(m)];
                if (metab > 0.01)
                {
                    SCFAFILE << "," << metab;
                }
                else
                {
                    SCFAFILE << "," << 0;
                }
            }
            SCFAFILE << std::endl;
        }
    }
    SCFAFILE.flush();
}

void Gut_Output::Write_BacteriaSummary(Grid &g)
{
    /*
     This method writes to the output file "BacteriaSummary.csv"
     with the prevalances of all bacteria present, scaled by volume
     Placed in the simulationdirectory directly for easy access
     */
    // map <string,int> bacteriamap;
    for (int j = 1; j <= g.ncol; j++)
    {
        for (int i = 0; i <= g.nrow + 1; i++)
        {
            string species = g.Cellfield[i][j][0].species;
            if (bacteriamap.find(species) == bacteriamap.end())
            {
                bacteriamap[species] = g.Cellfield[i][j][0].volume;
            }
            else
            {
                bacteriamap[species] += g.Cellfield[i][j][0].volume;
            }
        }
    }
    BACTERIAFILE << simultime;
    for (std::map<string, double>::iterator it = bacteriamap.begin(); it != bacteriamap.end(); ++it)
    {
        BACTERIAFILE << "," << it->first << "," << it->second;
    }
    BACTERIAFILE << std::endl;
    BACTERIAFILE.flush();

    // Set bacteria values back to 0
    for (const auto &p : bacteriamap)
    {
        bacteriamap[p.first] = 0;
    }
}

void Gut_Output::Write_Metabolites(Grid &g)
{
    /* Print total metabolite amounts above trivial
     */
    map<string, double> totmetabolitemap;
    for (int i = 1; i <= g.nrow; i++)
    {
        for (int j = 1; j <= g.ncol; j++)
        {
            for (auto const &x : g.exchange_metab_dict)
            {
                double metabconc = g.metfield[i][j][x.second];
                if (metabconc > 0)
                {
                    if ((x.first != "h2o"))
                    {
                        if (totmetabolitemap.find(x.first) == totmetabolitemap.end())
                        {
                            totmetabolitemap[x.first] = metabconc;
                        }
                        else
                        {
                            totmetabolitemap[x.first] += metabconc;
                        }
                    }
                }
            }
        }
    }
    // std::cout << "total number of limited metabolites " << totmetabolitemap.size();
    // std::cout << std::endl;
    std::vector<pair<double, std::string>> metabsvectorised;

    for (std::map<string, double>::iterator it = totmetabolitemap.begin(); it != totmetabolitemap.end(); ++it)
    {
        double metab = it->second;
        if (metab > 0.1)
        {
            metabsvectorised.push_back(make_pair(it->second, it->first));
        }
    }
    sort(metabsvectorised.rbegin(), metabsvectorised.rend());
    for (unsigned i = 0; i < metabsvectorised.size(); i++) // i will only be positive, and will be compared with another unsigned - the size
    {
        METABFILE << simultime << "," << metabsvectorised[i].first << "," << metabsvectorised[i].second;
        METABFILE << std::endl;
    }
    METABFILE.flush();
}

void Gut_Output::Write_Reactions(Grid &g)
{
    // Print the reactions used by each bacterium

    REACTFILE << "reactions at " << simultime << std::endl;
    for (int i = 1; i <= g.nrow; i++)
    {
        for (int j = 1; j <= g.ncol; j++)
        {
            if (g.Cellfield[i][j][0].species == "none")
            {
                continue;
            }
            int totreact = 0;
            double totflux = 0;
            // std::cout << g.Cellfield[i][j][0].species << " " << i << " " << j << " " << std::endl;
            REACTFILE << g.Cellfield[i][j][0].species << " " << i << " " << j << " " << g.Cellfield[i][j][0].growthrate * par.growthmod << std::endl;
            for (auto const &x : g.Cellfield[i][j][0].react_dict)
            {
                if (g.Cellfield[i][j][0].used[x.second - 1] == 1)
                {
                    if (g.Cellfield[i][j][0].flux[x.second - 1] > 0.0000000000000001)
                    { // 0.000001){
                        totreact += 1;
                        // std::cout << x.first << std::endl;
                        REACTFILE << x.first << "," << g.Cellfield[i][j][0].flux[x.second - 1] << std::endl;
                        totflux += g.Cellfield[i][j][0].flux[x.second - 1];
                    }
                }
            }

            // if(totflux>0){
            // std::cout<<"Total flux here is "<<totflux<<std::endl;
            std::cout << "Total number of reactions is " << totreact << " for " << g.Cellfield[i][j][0].species << std::endl;
            // std::cout<<"Total growth is "<<g.Cellfield[i][j][0].growthrate <<std::endl;
            //           }

            REACTFILE << std::endl;
        }
    }
    REACTFILE.flush();

    // An alternate way of recording reactions that can be analysed systematically, but is not human-readable
    /*
        for (int i = 1; i <= g.nrow; i++)
        {
            for (int j = 1; j <= g.ncol; j++)
            {
                if(g.Cellfield[i][j][0].volume>0){
                //std::cout << g.Cellfield[i][j][0].species << " " << i << " " << j << " " << std::endl;
                REACTFILE << simultime <<","<< g.Cellfield[i][j][0].species<<","<<g.Cellfield[i][j][0].volume << "," << g.Cellfield[i][j][0].growthrate;
                for (auto const &x : g.Cellfield[i][j][0].react_dict)
                {
                    if (g.Cellfield[i][j][0].used[x.second - 1] == 1)
                    {
                        if(g.Cellfield[i][j][0].flux[x.second - 1]>0.00001){
                            //std::cout << x.first << std::endl;
                            REACTFILE <<","<< x.second;
                        }
                    }
                }
                REACTFILE << std::endl;
            }
            }
        }
        REACTFILE.flush();
    */
}

void Gut_Output::Write_Energy(Grid &g, std::map<std::string, double> gibbsmap)
{
    // Get energy entering the system
    double inenergy = 0;
    for (auto const &x : g.inputmetmap)
    {
        inenergy += (gibbsmap[x.first] * x.second);
    }
    std::cout << "Energy that entered the system is " << inenergy << " mJ" << std::endl;

    // Get energy of everything currently in system
    double systemenergy = 0;
    for (int i = 1; i <= g.nrow; i++)
    {
        for (int j = 1; j <= g.ncol; j++)
        {
            for (auto const &x : g.exchange_metab_dict)
            {
                if (g.metfield[i][j][x.second] > 0.0000000001)
                {
                    systemenergy += (gibbsmap[x.first] * g.metfield[i][j][x.second]);
                }
            }
        }
    }
    std::cout << "Energy currently in system is " << systemenergy << " mJ" << std::endl;

    // Get system output
    double outenergy = 0;
    for (auto const &x : g.fecesmetmap)
    {
        outenergy += (gibbsmap[x.first] * x.second);
    }
    std::cout << "Energy that left the system is " << outenergy << " mJ" << std::endl;

    TOTENERGYFILE << simultime << "," << std::setprecision(std::numeric_limits<double>::digits10 + 1) << inenergy << "," << systemenergy << "," << outenergy << std::endl;
    TOTENERGYFILE.flush();
}

void Gut_Output::Write_Elements(Grid &g, std::map<std::string, double> carbonmap, std::map<std::string, double> hydrogenmap, std::map<std::string, double> nitrogenmap, std::map<std::string, double> oxygenmap)
{
    // Get energy entering the system
    double incarbon = 0;
    double inhydrogen = 0;
    double innitrogen = 0;
    double inoxygen = 0;
    for (auto const &x : g.inputmetmap)
    {
        incarbon += (carbonmap[x.first] * x.second);
        inhydrogen += (hydrogenmap[x.first] * x.second);
        innitrogen += (nitrogenmap[x.first] * x.second);
        inoxygen += oxygenmap[x.first] * x.second;
    }
    std::cout << "Carbon that entered the system is " << incarbon << " um" << std::endl;

    // Get carbon of everything currently in system
    double systemcarbon = 0;
    double systemhydrogen = 0;
    double systemnitrogen = 0;
    double systemoxygen = 0;
    for (int i = 1; i <= g.nrow; i++)
    {
        for (int j = 1; j <= g.ncol; j++)
        {
            for (auto const &x : g.exchange_metab_dict)
            {
                if (g.metfield[i][j][x.second] > 0.0000000001)
                {
                    systemcarbon += (carbonmap[x.first] * g.metfield[i][j][x.second]);
                    systemhydrogen += (hydrogenmap[x.first] * g.metfield[i][j][x.second]);
                    systemnitrogen += (nitrogenmap[x.first] * g.metfield[i][j][x.second]);
                    systemoxygen += (oxygenmap[x.first] * g.metfield[i][j][x.second]);
                }
            }
        }
    }
    std::cout << "Carbon currently in system is " << systemcarbon << " um" << std::endl;

    // Get system output
    double outcarbon = 0;
    double outhydrogen = 0;
    double outnitrogen = 0;
    double outoxygen = 0;
    for (auto const &x : g.fecesmetmap)
    {
        outcarbon += (carbonmap[x.first] * x.second);
        outhydrogen += (hydrogenmap[x.first] * x.second);
        outnitrogen += (nitrogenmap[x.first] * x.second);
        outoxygen += (oxygenmap[x.first] * x.second);
    }
    std::cout << "Carbon that left the system is " << outcarbon << " uM" << std::endl;
    TOTELEMENTFILE << simultime << "," << std::setprecision(std::numeric_limits<double>::digits10 + 1) << incarbon << "," << systemcarbon << "," << outcarbon << "," << inhydrogen << "," << systemhydrogen << "," << outhydrogen << "," << innitrogen << "," << systemnitrogen << "," << outnitrogen << "," << inoxygen << "," << systemoxygen << "," << outoxygen << std::endl;
    TOTELEMENTFILE.flush();
}

void Gut_Output::Write_Feces(Grid &g)
{
    for (auto const &x : g.fecesbacmap)
    {
        FECESBACFILE << simultime << "," << x.first << "," << x.second << std::endl;
    }
    g.fecesbacmap.clear();
    FECESBACFILE.flush();
    std::vector<pair<double, std::string>> fecesvectorised;

    for (auto const &x : g.fecesmetmap)
    {
        double metab = x.second;
        if (metab > 0.000001)
        {
            fecesvectorised.push_back(make_pair(x.second, x.first));
        }
    }

    sort(fecesvectorised.rbegin(), fecesvectorised.rend());
    for (unsigned i = 0; i < fecesvectorised.size(); i++) // n will only be positive, and will be compared with another unsigned - the size
    {
        FECESMETFILE << simultime << "," << fecesvectorised[i].first << "," << fecesvectorised[i].second;
        FECESMETFILE << std::endl;
    }
    FECESMETFILE.flush();
}

void Gut_Output::Write_Host_Uptake(Grid &g)
{
    map<string, double> readablemetabolitemap;
    for (auto const &x : g.exchange_metab_dict)
    {
        if (x.first != "h2o")
        {
            double metab = g.uptakemetmap[x.second];
            if (metab > 0.001)
            {
                UPTAKEMETFILE << simultime << "," << x.first << "," << metab << std::endl;
            }
        }
    }
    g.uptakemetmap.clear();
    UPTAKEMETFILE.flush();
}

void Gut_Output::Write_Params()
{

    PARAMFILE << "seed " << par.randseed << std::endl;
    PARAMFILE << "lactose input " << par.foodinlactose << std::endl;
    PARAMFILE << "oxygen input " << par.oxin << std::endl;
    PARAMFILE << "initial oxygen " << par.initialox << std::endl;
    PARAMFILE << "cell drift " << par.celldrift << std::endl;
    PARAMFILE << "time between food pulses in hours " << par.timefood << std::endl;
    PARAMFILE << "time steps per hour " << par.timescale << std::endl;
    PARAMFILE << "metabolite flow" << par.metdrift << std::endl;
    PARAMFILE << "placement chance (location depends on placement method) = " << par.placechance << std::endl;
    PARAMFILE << "Maximum metabolic uptake = " << par.max_metup << std::endl;
    PARAMFILE << "flux weight = " << par.fluxweight << std::endl;
    PARAMFILE << "Initial placement chance = " << par.INITDENS << std::endl;
    PARAMFILE << "Diffusion chance = " << par.diffconst << std::endl;
    PARAMFILE << "CHOSEN TRACKED METABOLITES" << std::endl;
    for (std::string m : par.outputmetabs)
    {
        PARAMFILE << m << std::endl;
    }
    PARAMFILE.flush();
}
