#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include "BactParam.h"
#include "Cell.h"
#include "Grid.h"
#include "Gut_Output.h"
#include "Graphics.h"
#include "parameter.h"
#include "project.h"

using namespace std;

extern Par par;

std::string simulationdir;
Gut_Output OUTPUT;

std::map<std::string, double> gibbsmap;
std::map<std::string, double> carbonmap;
std::map<std::string, double> hydrogenmap;
std::map<std::string, double> nitrogenmap;
std::map<std::string, double> oxygenmap;

int ancestorcount = 0;
int simultime;

int main(int argc, char *argv[])
{
    std::cout << "started" << std::endl;
    // The command line arguments are parsed here
    par.runid = atoi(argv[1]);
    if (par.runid > 100)
    {
        par.randseed = par.runid;
    }
    else
    {
        srand(time(0) + par.runid);
        par.randseed = rand() + 150 - par.runid;
    }
    std::cout << "Random seed is: " << par.randseed << std::endl;
    int experimentmode = 0; // Default, might be changed later
    int n_iterations = atoi(argv[2]);
    par.foodinlactate = atof(argv[3]);
    par.foodinGOS = atof(argv[4]);
    par.foodinfl = atof(argv[5]);
    par.foodinlactose = atof(argv[6]);
    par.oxin = atof(argv[7]);
    par.initialox = atof(argv[8]);
    par.metdrift = atof(argv[9]);
    par.ncol = atof(argv[10]);
    par.nrow = atof(argv[11]);
    par.biomassreplace = atof(argv[12]);
    if (atof(argv[13]) > 0)
    {
        par.alwayssaveenergy = true;
    }
    par.celldrift = atof(argv[14]);
    if (par.celldrift == 4)
    { // Special case: enable in vitro conditions
        par.frequentdata = true;
        par.shuffle = true;
        par.intakeincrease = -1000;
        par.metdrift = 0;
        par.celldrift = 0;
        par.ncol = 225;
        par.nrow = 9;
    }
    par.bacplacetype = atof(argv[15]);

    // Read in like this if we're using an extra argument to assign a filename
    if (par.bacplacetype == 3)
    {
        par.bacplacefile = argv[16];
        par.immigration = std::floor(atof(argv[17]));
        par.placechance = std::fmod(atof(argv[17]), 1.0); // fmod is % but for doubles, essentially
        par.deathsetting = atof(argv[18]);

        if (par.deathsetting == 2)
        {
            par.pdeath = 2;
            par.pdeath_base = 0.025;
        }
        else
        {
            par.pdeath = 0;
            par.pdeath_base = par.deathsetting;
        }
        par.growthmod = atof(argv[19]);
        par.bactmove = atof(argv[20]);
        par.initialvolume = atof(argv[21]);
        experimentmode = atoi(argv[22]);
        par.fluxweight = atof(argv[23]);
        par.diffconst = atof(argv[24]);
        par.cellfrac = atof(argv[25]);
        par.testname = argv[26];
    }
    // Read in like this if we're not using an extra argument to assign a filename
    else
    {
        par.immigration = std::floor(atof(argv[16]));
        par.placechance = std::fmod(atof(argv[16]), 1.0); // fmod is % but for doubles, essentially
        par.deathsetting = atof(argv[17]);

        if (par.deathsetting == 2)
        {
            par.pdeath = 2;
            par.pdeath_base = 0.025;
        }
        else
        {
            par.pdeath = 0;
            par.pdeath_base = par.deathsetting;
        }
        par.growthmod = atof(argv[18]);
        par.bactmove = atof(argv[19]);
        par.initialvolume = atof(argv[20]);
        experimentmode = atoi(argv[21]);
        par.fluxweight = atof(argv[22]);
        par.diffconst = atof(argv[23]);
        par.cellfrac = atof(argv[24]);
        par.testname = argv[25];
    }
    // Set our diffussion correctly (more steps if a high value has been entered)
    int diffusesteps = 1;
    if (par.diffconst > 4)
    {
        diffusesteps = par.diffconst / 4;
        par.diffconst = 2.85;
    }
    //Do the same for the bacterial diffusion (through swaps)
    int swapsteps = 1;  
    if(par.bactmove>1){
    swapsteps = par.bactmove;
    }
    if(par.bactmove<1){
        par.maxswaps = 740*par.bactmove;
    }

    // Assign the right experiment mode to knock out certain things for specific experiments
    switch (experimentmode)
    {
    case 0:
        break;
    case 1:
        par.knockoutlactate = true;
        break;
    case 2:
        par.knockoutxfp = true;
        break;
    case 3:
        par.knockoutnonbiflactose = true;
        break;
    case 4:
        par.knockoutnonbiflactoselate = true;
        break;
    case 5:
        par.knockoutenterooxygen = true;
        break;
    case 6:
        par.knockoutbiflactate = true;
        break;
    case 7:
        par.knockoutecollactate = true;
        break;
    case 8:
        par.knockoutbifandecollactate = true;
        break;
    case 9:
        par.diminishoxygenrelease = true;
        break;
    case 10:
        swapsteps = 10;
        break;
    case 11:
        par.maxswaps = 74;
        break;
    case 12:
        swapsteps = 5;
        break;
    case 13:
        par.maxswaps = 148;
        break;
    case 14:
        par.knockoutbiflactateup = true;
        break;
    case 15:
        par.knockoutbutyro12ppd = true;
        break;
    case 16:
        par.knockout12ppd = true;
        break;
    case 17:
        par.knockoutbutyrolcts = true;
        break;
    case 18:
        par.knockoutbutyrolactate = true;
        break;
     case 19:
        par.knockoutbutyrolactoselate= true;
        break;
    }
    std::cout << "Knockouts enabled: " << par.knockoutlactate << par.knockoutxfp << par.knockoutnonbiflactose << par.knockoutnonbiflactoselate << par.knockoutenterooxygen << par.knockoutbiflactate << par.knockoutecollactate << par.knockoutbifandecollactate << std::endl;

    simulationdir = par.testname;
    std::string commandLineStr = "";
    for (int i = 1; i < argc - 1; i++)
    {
        simulationdir += "_";
        simulationdir += argv[i];
    }

    // Make a map of all available metabolic gibbs values
    // This is used for some sanity checking on reactions.
    std::ifstream rawgibbs;
    rawgibbs.open("energiesdatafile.csv");
    std::string nextline;
    if (rawgibbs.is_open())
    {
        while (getline(rawgibbs, nextline))
        {
            // get id
            int index = nextline.find(";");
            std::string id = nextline.substr(0, index);
            std::string gibbs = nextline.substr(index + 1, nextline.length());

            // add value to map if known and not already in map
            if (gibbs != "UNKNOWN")
            {
                if (gibbsmap.find(id) == gibbsmap.end())
                {
                    gibbsmap[id] = stod(gibbs);
                }
            }
        }
        rawgibbs.close();
    }

    // Fill b_params with the different bacteria.
    QVector<BactParam> b_params;
    std::cout << "Initializing species" << std::endl;

   //Species based on Backhed Newborn list, arbitrary order
    par.knockoutbifandecollactate = true;
    AGORAParam bp1 ("MODEL_Bvul.xml",1,par.fluxweight,gibbsmap);
    AGORAParam bp2 ("MODEL_Sora.xml", 2, par.fluxweight,gibbsmap);
    AGORAParam bp3 ("MODEL_BlongInf.xml", 3, par.fluxweight,gibbsmap);
    AGORAParam bp4 ("MODEL_Sepi.xml", 4, par.fluxweight,gibbsmap); 
    AGORAParam bp5 ("MODEL_Cbut.xml", 5, par.fluxweight,gibbsmap); //butyrate
    AGORAParam bp6 ("MODEL_Gmor.xml", 6, par.fluxweight,gibbsmap); 
    AGORAParam bp7 ("MODEL_Rmuc.xml", 7, par.fluxweight,gibbsmap);
    AGORAParam bp8 ("MODEL_Ecol.xml", 8, par.fluxweight,gibbsmap);
    AGORAParam bp9 ("MODEL_Caer.xml", 9, par.fluxweight,gibbsmap);
    AGORAParam bp10 ("MODEL_Efae.xml", 10, par.fluxweight,gibbsmap);
    AGORAParam bp11 ("MODEL_Lgas.xml", 11, par.fluxweight,gibbsmap);
    AGORAParam bp12 ("MODEL_Pdis.xml", 12, par.fluxweight,gibbsmap);
    AGORAParam bp13 ("MODEL_Rgna.xml", 13, par.fluxweight,gibbsmap);
    AGORAParam bp14 ("MODEL_Pacn.xml",14,par.fluxweight,gibbsmap); 
    AGORAParam bp15 ("MODEL_Hpar.xml", 15, par.fluxweight,gibbsmap); 
    AGORAParam bp16 ("MODEL_Ehal.xml", 16, par.fluxweight,gibbsmap); //butyrate
    AGORAParam bp17 ("MODEL_Vdis.xml", 17, par.fluxweight,gibbsmap);
    AGORAParam bp18 ("MODEL_EYY.xml", 18, par.fluxweight,gibbsmap); 
    AGORAParam bp19 ("MODEL_Rinu.xml", 19, par.fluxweight,gibbsmap); //butyrate
    AGORAParam bp20 ("MODEL_BlongLong.xml", 20, par.fluxweight,gibbsmap);

    std::cout << "initializing parameters" << std::endl;
    b_params.append(bp1);
    b_params.append(bp2);
    b_params.append(bp3);
    b_params.append(bp4);
    b_params.append(bp5);
    b_params.append(bp6);
    b_params.append(bp7);
    b_params.append(bp8);
    b_params.append(bp9);
    b_params.append(bp10);
    b_params.append(bp11);
    b_params.append(bp12);
    b_params.append(bp13);
    b_params.append(bp14);
    b_params.append(bp15);
    b_params.append(bp16);
    b_params.append(bp17);
    b_params.append(bp18);
    b_params.append(bp19);
    b_params.append(bp20);
    
    bp3.remove_reaction("R_LDH_L");
    bp8.remove_reaction("R_LDH_L");
    bp20.remove_reaction("R_LDH_L");
    std::cout << "parameters done" << std::endl;

    // Make a vector containing all the metabolites listed in the SBML models (union, not intersection).
    QVector<std::string> e_metab;
    Exchange exch;
    for (auto const &bp : b_params)
    {
        for (auto const &x : bp.exchange_dict)
        {
            if (nitrogenmap.find(x.first) == nitrogenmap.end())
            {
                exch = x.second;
                nitrogenmap[x.first] = exch.nitrogen;
            }
            if (carbonmap.find(x.first) == carbonmap.end())
            {
                exch = x.second;
                carbonmap[x.first] = exch.carbon;
            }
            if (hydrogenmap.find(x.first) == hydrogenmap.end())
            {
                exch = x.second;
                hydrogenmap[x.first] = exch.hydrogen;
            }
            if (oxygenmap.find(x.first) == oxygenmap.end())
            {
                exch = x.second;
                oxygenmap[x.first] = exch.oxygen;
            }
            if (not e_metab.contains(x.first))
                e_metab.append(x.first);
        }
    }
    std::cout << "vector done" << std::endl;

    // workaround because mkdir is actually deprecated
    char simdir[simulationdir.size() + 1];
    strcpy(simdir, simulationdir.c_str());
    if (mkdir(simdir, 0777) == -1) // creating a directory
    {
        cerr << "Directory " << simulationdir << " already exists" << endl; // This is also fine
    }
    std::cout << "Creating grid" << std::endl;
    Grid My_Grid(par.nrow, par.ncol, e_metab); // Creating the grid with specified number of columns and rows and all metabolites in e_metab.
    std::cout << "Filling grid" << std::endl;
    My_Grid.Init(b_params); // Filling it with bacteria from b_params
    if (par.bacplacetype == 3)
    {
        My_Grid.CellShuffle(); // Shuffle cells if placed from file
    }
    std::cout << "Creating output folder" << std::endl;
    OUTPUT.Init(); // Create data folder with output files
    std::cout << "Initializing graphics" << std::endl;
    // Create folder for  movie and prepare to generate pngs
    Graphics my_graphics(My_Grid);
    my_graphics.Init_Movie();
    OUTPUT.Write_Params();

    time_t t = time(0); // get time now
    struct tm *now = localtime(&t);
    std::cout << "Start of simulation: " << (now->tm_hour) << ':' << (now->tm_min) << ':' << (now->tm_sec) << std::endl;
    simultime = 0;
    OUTPUT.Write_BacteriaSummary(My_Grid);
    OUTPUT.Write_Output(My_Grid);
    OUTPUT.Write_SCFA(My_Grid);
    OUTPUT.Write_Metabolites(My_Grid);
    for (simultime = 1; simultime <= n_iterations; simultime++) // 1 h = 20 timesteps
    {
        std::cout << "step " << simultime << " of " << simulationdir << std::endl;
        if ((simultime % (int)(par.timefood * par.timescale) == 1))
        {
            My_Grid.MixedSugarInflux();
        }
        My_Grid.HostInflux();               // Add things at the gut wall. Under normal parameters, nothing is added here.
        My_Grid.Set_Unlimited_Metabolite(); // determine which metabolites are unlimited, set those concentrations to 9999

        if (((simultime % (int)(par.timefood * par.timescale) == 1)) || par.frequentdata)
        {
            // Do these recordings at every feeding
            OUTPUT.Write_Energy(My_Grid, gibbsmap);
            OUTPUT.Write_Feces(My_Grid);
            OUTPUT.Write_Elements(My_Grid, carbonmap, hydrogenmap, nitrogenmap, oxygenmap);
            // Clear the relevant lists after they have been recorded
            My_Grid.fecesmetmap.clear();
            My_Grid.inputmetmap.clear();
            // OUTPUT.Write_Host_Uptake(My_Grid); #Enable if doing host uptake
        }

        // Make conditions well-mixed, if that is turned on
        if (par.shuffle)
        {
            My_Grid.CellShuffle();
            My_Grid.Mix_Metabolites();
        }

        // Do Metabolism
        //std::cout << "metabolism " << simultime << std::endl;
        My_Grid.Step_Metabolism(b_params);
        // growth, uptake and excretion rates of all Cells have been calculated: now Cell division, Cell death and metabolite diffusion are done
        if (par.metdrift < 1000 && par.metdrift > 0)
        {
            if (par.metdrift > 10)
            {
                My_Grid.DriftCells(b_params); // Cells and metabolites all drift to the right
                if (par.metdrift > 20)
                {                                 // If our drift is very high, move again
                    My_Grid.DriftCells(b_params); // Cells and metabolites all drift to the right
                }
            }
            else
            {
                if (simultime % (int)(par.timescale / par.metdrift) == 2)
                {
                    My_Grid.DriftCells(b_params); // Cells and metabolites all drift to the right
                }
            }
        }
        //std::cout << "celldiffusion " << simultime << std::endl;
        for (int i = 0; i < swapsteps; i++)
        {
            My_Grid.Diffuse_Cells(); //Kawasaki-Ising dynamics
            std::cout << "Bacterial diffusion  " << simultime << std::endl;
        }

        for (int i = 0; i < diffusesteps; i++)
        {
            My_Grid.Diffuse_Metabolites(par.diffconst / par.timescale);
            std::cout << "Diffusion  " << simultime << std::endl;
        }
        //std::cout << "placement " << simultime << std::endl;
        My_Grid.Place_Bacterium(b_params);
        //std::cout << "reproduction " << simultime << std::endl;
        My_Grid.Grid_Die_Reproduce();
        //std::cout << "data output " << simultime << std::endl;
        // append data to output files
        OUTPUT.Write_Metabolites(My_Grid);
        OUTPUT.Write_BacteriaSummary(My_Grid);

        OUTPUT.Write_SCFA(My_Grid);
        // This is for the exchange.dat file, recording relevant fluxes. This should be done every step if you want good data, but that does make for a very large file
        OUTPUT.Write_Output(My_Grid);
        // Make png of current iteration every step
        my_graphics.Create_Planes(My_Grid);
        my_graphics.Write_Movie(My_Grid);

        // Write reactions used per population to file every 400 steps
        if (simultime % 400 == 1)
        {
            OUTPUT.Write_Reactions(My_Grid);
        }

    } // end of time-loop

    time_t t2 = time(0); // get current time
    struct tm *now2 = localtime(&t2);
    std::cout << "End of simulation: " << (now2->tm_hour) << ':' << (now2->tm_min) << ':' << (now2->tm_sec) << std::endl;
    std::cout << "GLPK will now abort. This is normal behaviour." << std::endl;

    return 0;
}
