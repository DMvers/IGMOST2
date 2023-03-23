#include <sys/stat.h>
#include <sys/types.h>
#include <pngwriter.h>

#include <sys/stat.h>
#include <sys/types.h>
#include "parameter.h"
#include "Graphics.h"

/* The Graphics class produces images that visualise the grid (the bacteria and
 * the concentrations of different metabolites), these images can be turned into a video.
 */

// extern int nrow,ncol,scale;
extern Par par;
extern std::string simulationdir;

Graphics::Graphics(const Grid &g)
{
    /* This method constructs the Graphics class on the basis of an instance of the Grid. */
    nrow = g.nrow;
    ncol = g.ncol;
    numfields = 45;
    fields_display = new int **[numfields];
    for (int k = 0; k < numfields; ++k)
    {
        fields_display[k] = new int *[nrow + 2];
        for (int i = 0; i < nrow + 2; ++i)
        {
            fields_display[k][i] = new int[ncol + 2];
            for (int j = 0; j < ncol + 2; j++)
            {
                fields_display[k][i][j] = 0;
            }
        }
    }
}

Graphics::Graphics(const Graphics &g)
{
    /* This method constructs the Graphics class on the basis of a previous instance of Graphics. */
    numfields = g.numfields;
    countmovie = g.countmovie;
    fields_display = new int **[numfields];
    for (int k = 0; k < numfields; ++k)
    {
        fields_display[k] = new int *[nrow + 2];
        for (int i = 0; i < nrow + 2; ++i)
        {
            fields_display[k][i] = new int[ncol + 2];
            for (int j = 0; j < ncol + 2; j++)
            {
                fields_display[k][i][j] = g.fields_display[k][i][j];
                ;
            }
        }
    }
}

Graphics::~Graphics()
{
    /* The destructor of the Graphics class. */
    for (int k = 0; k < numfields; k++)
    {
        for (int i = 0; i < nrow + 2; i++)
            delete[] fields_display[k][i];
        delete[] fields_display[k];
    }
    delete[] fields_display;
}

Graphics &Graphics::operator=(const Graphics &g)
{
    /* Making a reference of the Graphics class.  */
    if (this == &g)
    {
        return *this;
    }

    for (int k = 0; k < numfields; k++)
    {
        for (int i = 0; i < nrow + 2; i++)
            delete[] fields_display[k][i];
        delete[] fields_display[k];
    }

    fields_display = new int **[numfields];
    for (int k = 0; k < numfields; ++k)
    {
        fields_display[k] = new int *[nrow + 2];
        for (int i = 0; i < nrow + 2; ++i)
        {
            fields_display[k][i] = new int[ncol + 2];
            for (int j = 0; j < ncol + 2; j++)
            {
                fields_display[k][i][j] = g.fields_display[k][i][j];
            }
        }
    }
    numfields = g.numfields;
    countmovie = g.countmovie;
    return *this;
}

void Graphics::Init_Movie()
{
    /* Initialisation of the Graphics class. The movie directory is created. */
    moviedir = simulationdir;
    moviedir += par.moviedir;
    // workaround because mkdir is actually deprecated
    char mdir[moviedir.size() + 1];
    strcpy(mdir, moviedir.c_str());
    if (mkdir(mdir, 0777) == -1) // creating a directory
    {
        cerr << "Directory " << moviedir << " already exists" << endl;
    }

    countmovie = 0;
}

void Graphics::Create_Planes(const Grid &g)
{
    /* This method creates the different visualisations of the grid for the movie.
     * Such as the bacteria and rows for the concentrations of metabolites
     */
    for (int i = 1; i <= nrow; i++)
    {
        for (int j = 1; j <= ncol; j++)
        {
            Celllayer = new Cell[par.maxcells];
            for (int k = 0; k < par.maxcells; k++)
            {
                Celllayer[k].exist = false;
            }
            for (int k = 0; k < par.maxcells; k++)
            {
                if (g.Cellfield[i][j][k].exist)
                {
                    Celllayer[k] = g.Cellfield[i][j][k];
                }
            }
            std::map<std::string, int> emd = g.exchange_metab_dict;
            int smallmagnification = 220;

            // fields_display 1:20 reserved for metabolites
            fields_display[1][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("12ppd_S")])); // 12ppd_S
            fields_display[2][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("lcts")]));  // lactose
            fields_display[3][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("ac")]));    // acetate
            fields_display[4][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("but")]));    //butyrate
            fields_display[6][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("h2")]));    // hydrogen
            fields_display[7][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("etoh")]));  // ethanol
            fields_display[8][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("succ")]));  // succinate
            fields_display[9][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("2FuLa")]));    // 2'-FL
            fields_display[10][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("gal")]));  // galactose
            fields_display[11][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("co2")]));  // co2
            fields_display[12][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("GOS3")]));    // GOS DP3
            fields_display[13][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("GOS4")]));    // GOS DP4
            fields_display[14][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("GOS5")]));    // GOS DP5
            fields_display[15][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("for")]));  // formate
            fields_display[16][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("lac_L")] +
                                                          smallmagnification * g.metfield[i][j][emd.at("lac_D")])); // l-lactate + d-lactate
            fields_display[17][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("for")]));   // formate
            fields_display[18][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("o2")]));    // oxygen
            fields_display[19][i][j] = min(206, 3 + (int)(smallmagnification * g.metfield[i][j][emd.at("ppa")]));   // propionate

            int k = 0;
            for (int n = 21; n < 41; n++)
            {
                fields_display[n][i][j] = 0;
            }
            fields_display[20 + Celllayer[k].indexnr][i][j] = max(30, 5 + (int)(600 * Celllayer[k].growthrate));

            delete[] Celllayer;
        }
    }
}

void Graphics::Write_Movie(const Grid &g)
{
    /* This method renders the actual PNGs for the movie. */
    char filename[200];
    char mdir[moviedir.size() + 1];
    strcpy(mdir, moviedir.c_str());
    sprintf(filename, "%s/%.5d.png", mdir, countmovie);
    pngwriter image(ncol, 14 * nrow + 14, 1.0, filename); // Number before nrow should be 1 higher than the highest one below
    for (int i = 1; i <= nrow; i++)
    {
        for (int j = 1; j <= ncol; j++)
        {
              int whitespecies = fields_display[22][i][j] + fields_display[24][i][j] + fields_display[26][i][j] +
                               fields_display[27][i][j] + fields_display[29][i][j] + fields_display[30][i][j] + fields_display[31][i][j] + fields_display[32][i][j] + fields_display[33][i][j]
                               + fields_display[34][i][j] + fields_display[35][i][j]+ fields_display[37][i][j]+ fields_display[38][i][j]; 

            int redspecies = fields_display[23][i][j] + fields_display[40][i][j];  // Bifidobacteria
            int greenspecies = fields_display[25][i][j] + fields_display[36][i][j] + fields_display[39][i][j];   // Butyrogenics
            int orangespecies = 0;      // unused
            int bluespecies = fields_display[28][i][j];  // E. coli
            int pinkspecies = fields_display[21][i][j];  // Bacteroides vulgatus

            if (whitespecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, min(1300 * whitespecies, 65534), min(1300 * whitespecies, 65534), min(1300 * whitespecies, 65534));
            }
            else if (redspecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, min(1300 * redspecies, 65534), 0, 0);
            }
            else if (bluespecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, 0, 0, min(1300 * bluespecies, 65534));
            }
            else if (pinkspecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, min(1300 * pinkspecies, 65534), 0, min(1300 * pinkspecies, 65534));
            }
            else if (greenspecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, 0, min(1300 * greenspecies, 65534), 0);
            }
            else if (orangespecies > 0)
            {
                image.plot(j, 13 * nrow + 13 + i, min(1300 * orangespecies, 65534), min(1300 * orangespecies, 65534) / 2, 0);
            }
            else
            {
                image.plot(j, 13 * nrow + 13 + i, 0, 0, 0); // No bacteria here
            }

            // Settings for exploring lactose and oxygen
            image.plot(j, 12 * nrow + 12 + i, min(250 * fields_display[2][i][j], 65534), 0, 0);                                                                                    // lactose
            //image.plot(j, 11 * nrow + 11 + i, min(250 * fields_display[18][i][j], 65534), min(250 * fields_display[18][i][j], 65534), min(250 * fields_display[18][i][j], 65534)); // o2
            image.plot(j, 11 * nrow + 11 + i, min(250 * fields_display[9][i][j], 65534), 0, 0); // 2fl
            image.plot(j, 10 * nrow + 10 + i, 0, min(250 * fields_display[16][i][j], 65534), 0);                                                                                   // lactate (L and D combined)
            image.plot(j, 9 * nrow + 9 + i, 0, 0, min(250 * fields_display[3][i][j], 65534));                                                                                      // acetate
            image.plot(j, 8 * nrow + 8 + i, 0, 0, min(250 * fields_display[1][i][j], 65534));                                                                                      // ethanol
            image.plot(j, 7 * nrow + 7 + i, 0, 0, min(250 * fields_display[4][i][j], 65534));                                                                                     // formate
            image.plot(j, 6 * nrow + 6 + i, 0, 0, min(250 * fields_display[8][i][j], 65534));                                                                                      // Succinate
            image.plot(j, 5 * nrow + 5 + i, 0, min(250 * fields_display[11][i][j], 65534), min(250 * fields_display[11][i][j], 65534));                                            // co2
            image.plot(j, 4 * nrow + 4 + i, 0, min(250 * fields_display[6][i][j], 65534), min(250 * fields_display[6][i][j], 65534));                                              // hydrogen
            image.plot(j, 3 * nrow + 3 + i, 0, 0, min(250 * fields_display[19][i][j], 65534));                                                                                     // propionate
            image.plot(j, 2 * nrow + 2 + i, min(250 * fields_display[12][i][j], 65534), min(250 * fields_display[13][i][j], 65534), min(250 * fields_display[14][i][j], 65534));                                                                                                                              // empty
            image.plot(j, nrow + 1 + i, min(250 * fields_display[9][i][j], 65534), 0, 0);                                                                                                                                  // empty
            image.plot(j, i, 0, 0, 0);                                                                                                                                             // empty
        }
    }
    image.close();
    countmovie++;
}