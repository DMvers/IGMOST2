#include "project.h"
#include "Grid.h"
#include "parameter.h"

extern int simultime;
extern Par par;
extern int ancestorcount;
extern std::string simulationdir;

/** This class implements the grid which contains the bacteria and metabolite concentrations.
 */

Grid::Grid(int nr, int nc, QVector<std::string> exchange_metab)
{
	/* This method constructs the grid with dimension nr X nc.
	 * The vector exchange_metab contains the names of all the metabolites in the system (taken from the SBML models).
	 * A Cellfield is created with the bacteria.
	 * metfield contains the concentrations of metabolites, the metabolites are numbered.
	 * enzymefield contains a map of string to double for each grid position, containing a linkage (like galactose-glucose) and a strength (0-1)
	 */
	nrow = nr;
	ncol = nc;
	numexchangemet = exchange_metab.size();
	for (int i = 0; i < numexchangemet; ++i)
		exchange_metab_dict[exchange_metab.at(i)] = i;

	Cellfield = new Cell **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		Cellfield[i] = new Cell *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			Cellfield[i][j] = new Cell[par.maxcells];
		}
	}

	metfield = new double **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		metfield[i] = new double *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			metfield[i][j] = new double[numexchangemet];
			for (int k = 0; k < numexchangemet; k++)
			{
				metfield[i][j][k] = 0; // k specifies the metabolite
			}
		}
	}

	sumconc = new double[numexchangemet];
	for (int k = 0; k < numexchangemet; k++)
	{
		sumconc[k] = 0; // This array contains the average concentrations of the metabolites through the whole system.
	}
}

Grid::~Grid()
{
	/* Destructor of the Grid class. */
	for (int i = 0; i < nrow + 2; i++)
	{
		for (int j = 0; j < ncol + 2; j++)
		{
			delete[] metfield[i][j];
			delete[] Cellfield[i][j];
		}
		delete[] metfield[i];
		delete[] Cellfield[i];
	}
	delete[] Cellfield;
	delete[] metfield;

	delete[] sumconc;
}

Grid &Grid::operator=(const Grid &g)
{
	/* Reference operator for the Grid class. */
	if (this == &g)
	{
		return *this;
	}

	for (int i = 0; i < nrow + 2; i++)
	{
		for (int j = 0; j < ncol + 2; j++)
		{
			delete[] metfield[i][j];
			delete[] Cellfield[i][j];
		}
		delete[] metfield[i];
		delete[] Cellfield[i];
	}
	delete[] Cellfield;
	delete[] metfield;

	delete[] sumconc;

	nrow = g.nrow;
	ncol = g.ncol;
	numexchangemet = g.numexchangemet;
	exchange_metab_dict = g.exchange_metab_dict;

	Cellfield = new Cell **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		Cellfield[i] = new Cell *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			Cellfield[i][j] = new Cell[par.maxcells];
		}
	}

	metfield = new double **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		metfield[i] = new double *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			metfield[i][j] = new double[numexchangemet];
			for (int k = 0; k < numexchangemet; k++)
			{
				metfield[i][j][k] = 0;
			}
		}
	}

	sumconc = new double[numexchangemet];
	for (int k = 0; k < numexchangemet; k++)
	{
		sumconc[k] = 0;
	}

	return *this;
}

void Grid::Init(QVector<BactParam> bps)
{
	/* initializer for the Grid class.
	 * This method fills the Cellfield with bacteria according to a certain initial density and
	 * with an equal distribution of the different species.
	 */

	// Load bacplace.csv here if we're using that placement method
	std::vector<std::tuple<std::string, int>> bactoplace;
	unsigned currentbac = 0; // Can only be positive, and will be compared with a size
	if (par.bacplacetype == 3)
	{
		std::ifstream bacplacestream;
		bacplacestream.open(par.bacplacefile);
		std::string nextline;
		if (bacplacestream.is_open())
		{
			while (getline(bacplacestream, nextline))
			{
				// get name
				int index = nextline.find(";");
				std::string speciesname = nextline.substr(0, index);
				bacspecies.push_back(speciesname);
				// get count
				std::reverse(nextline.begin(), nextline.end());
				index = nextline.find(";");
				std::string baccount = nextline.substr(0, index);
				std::reverse(baccount.begin(), baccount.end());
				bactoplace.push_back(std::tuple<std::string, int>(speciesname, stoi(baccount)));
			}
		}
	}

	for (int i = 0; i <= nrow + 1; i++) // initialize cell field
	{
		for (int j = 0; j <= ncol + 1; j++)
		{
			for (int k = 0; k < par.maxcells; k++)
			{
				Cellfield[i][j][k] = Cell();
				if ((Uniform() < par.INITDENS || par.bacplacetype == 3) && j > 0 && j < ncol + 1 && i > 0 && i < nrow + 1 && par.bacplacetype != 0)
				{
					BactParam bp;
					ancestorcount++;
					if (par.bacplacetype == 1)
					{
						bp = bps.at(0); // Place only first bacterium
					}
					else if (par.bacplacetype == 2)
					{
						bp = bps.at(RandNum(bps.size()) - 1); // Randomly choose which species to create on this position.
					}
					if (par.bacplacetype == 3)
					{
						if (std::get<1>(bactoplace[currentbac]) == 0)
						{
							// std::cout << "currentbac " << currentbac << " of " << bactoplace.size() << std::endl;
							if (currentbac == bactoplace.size() - 1)
							{
								Cellfield[i][j][k] = Cell();
								// std::cout << "placed " << "nothing at " << i << " " << j << " " << k << std::endl;
								continue;
							}
							else
							{
								currentbac++;
							}
						}
						for (int b = 0; b < bps.length(); b++)
						{
							if (bps.at(b).speciesName == std::get<0>(bactoplace[currentbac]))
							{
								bp = bps.at(b); // Place only select bacterium
								Cellfield[i][j][k] = Cell(bp);
								Cellfield[i][j][k].Init(i, j);
								std::get<1>(bactoplace[currentbac])--;
								break;
							}
							if (b == bps.length() - 1)
							{
								std::cout << "Failed to find " << std::get<0>(bactoplace[currentbac]) << std::endl;
							}
						}
					}
					else
					{
						Cellfield[i][j][k] = Cell(bp);
						Cellfield[i][j][k].Init(i, j);
					}
				}
				else
				{
					Cellfield[i][j][k] = Cell(); // Create empty cell, no bacteria on this position.
				}
			}

			// Place initial oxygen if enabled
			if ((i != 0) && (i != nrow + 1))
			{
				if ((j != 0) && (j != ncol + 1))
				{
					if (inputmetmap.find("o2") == inputmetmap.end())
					{
						inputmetmap["o2"] = par.initialox;
					}
					else
					{
						inputmetmap["o2"] += par.initialox;
					}
					metfield[i][j][exchange_metab_dict.at("o2")] += par.initialox;
				}
			}
		}
	}
}

double Grid::CountCellNeighbour(int i, int j) const
{
	/* This method counts how many living neighbours (Moore neighbourhood) a bacteria at a certain position has. */
	// Returns a double 0<f<1 with how many of the possible neighbour positions are taken
	double count = 0;
	if (i == 1 && j == 1)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i + 1][j][k].exist + Cellfield[i][j + 1][k].exist + Cellfield[i + 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (4.0 * par.maxcells - 1);
	}
	else if (i == 1 && j == ncol)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i + 1][j][k].exist + Cellfield[i][j - 1][k].exist + Cellfield[i + 1][j - 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (4.0 * par.maxcells - 1);
	}
	else if (i == nrow && j == 1)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i][j + 1][k].exist + Cellfield[i - 1][j][k].exist + Cellfield[i - 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (4.0 * par.maxcells - 1);
	}
	else if (i == nrow && j == ncol)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i][j - 1][k].exist + Cellfield[i - 1][j][k].exist + Cellfield[i - 1][j - 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (4.0 * par.maxcells - 1);
	}
	else if (i == 1)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i][j - 1][k].exist + Cellfield[i][j + 1][k].exist + Cellfield[i + 1][j - 1][k].exist + Cellfield[i + 1][j][k].exist + Cellfield[i + 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (6.0 * par.maxcells - 1);
	}
	else if (i == nrow)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i][j - 1][k].exist + Cellfield[i][j + 1][k].exist + Cellfield[i - 1][j - 1][k].exist + Cellfield[i - 1][j][k].exist + Cellfield[i - 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (6.0 * par.maxcells - 1);
	}
	else if (j == 1)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i - 1][j][k].exist + Cellfield[i + 1][j][k].exist + Cellfield[i - 1][j + 1][k].exist + Cellfield[i][j + 1][k].exist + Cellfield[i + 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (6.0 * par.maxcells - 1);
	}
	else if (j == ncol)
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i - 1][j][k].exist + Cellfield[i + 1][j][k].exist + Cellfield[i - 1][j - 1][k].exist + Cellfield[i][j - 1][k].exist + Cellfield[i + 1][j - 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (6.0 * par.maxcells - 1);
	}
	else
	{
		for (int k = 0; k < par.maxcells; k++)
		{
			count += Cellfield[i - 1][j - 1][k].exist + Cellfield[i - 1][j][k].exist + Cellfield[i - 1][j + 1][k].exist + Cellfield[i][j - 1][k].exist + Cellfield[i][j + 1][k].exist + Cellfield[i + 1][j - 1][k].exist + Cellfield[i + 1][j][k].exist + Cellfield[i + 1][j + 1][k].exist + Cellfield[i][j][k].exist;
		}
		return (count - 1) / (9.0 * par.maxcells - 1); // added -1 to second value for consistency
	}
}

double *Grid::MetFieldNeighbour(int i, int j, int k, int n) const
{
	/* Get concentration of a certain metabolite (k) on a position neighbouring position (i, j) in a certain direction (n). */
	if (n == 1)
	{
		return &metfield[i][j + 1][k];
	}
	else if (n == 2)
	{
		return &metfield[i + 1][j][k];
	}
	else if (n == 3)
	{
		return &metfield[i][j - 1][k];
	}
	else if (n == 4)
	{
		return &metfield[i - 1][j][k];
	}
	return &metfield[i][j][k]; // Unused
}

void Grid::Update_Sumconc()
{
	/* This method calculates the concentrations of the different metabolites averaged over the entire grid. */
	for (int k = 0; k < numexchangemet; k++)
	{
		sumconc[k] = 0;
	}

	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= ncol; j++)
		{
			for (int k = 0; k < numexchangemet; k++)
			{
				sumconc[k] += metfield[i][j][k] / ((double)(nrow * ncol));
			}
		}
	}
}

void Grid::Diffuse_Cells()
{
	/* This method regulates the movement of the bacteria. */
	int swaps = 0;
	bool ***hasmoved;
	hasmoved = new bool **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		hasmoved[i] = new bool *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			hasmoved[i][j] = new bool[par.maxcells];
			for (int k = 0; k < par.maxcells; k++)
			{
				hasmoved[i][j][k] = 0;
			}
		}
	}
	int *randomgrid = new int[nrow * ncol];
	int *randomcell = new int[par.maxcells];
	for (int i = 0; i < nrow * ncol; i++)
	{
		randomgrid[i] = i;
	}
	Shufflevector(randomgrid, nrow * ncol);
	for (int tel = 0; tel <= nrow * ncol - 1; tel++)
	{
		int i = randomgrid[tel] / (ncol) + 1;
		int j = randomgrid[tel] % (ncol) + 1;
		int idir = 0;
		int jdir = 0;
		int k = 0;
		// determine direction
		int kdir = 0;

		int randomdirection = RandNum(8); //8 neighbours
		if(randomdirection<3){
			idir=-1;
		}
		if(randomdirection>4){
			idir=1;
		}
		if(randomdirection==0||randomdirection==3||randomdirection==5){
			jdir = -1;
		}
		if(randomdirection==2||randomdirection==4||randomdirection==7){
			jdir = 1;
		}
		//idir = RandNum(3) - 2; //+1,0, or -1
		//jdir = RandNum(3) - 2;

		if ((hasmoved[i][j][k] == false && hasmoved[i + idir][j + jdir][kdir] == false))
		{
			if (par.bactmove>0)
			{
				// std::cout << "moved"<<std::endl;

				if (i + idir != 0 && i + idir != nrow + 1) 
				{
					hasmoved[i][j][k] = true;
					hasmoved[i + idir][j + jdir][kdir] = true;
					if (j + jdir == ncol + 1) // moving out of the gut is only allowed at the rectal end
					{
						string species = Cellfield[i][j][0].species;

						if (fecesbacmap.find(species) == fecesbacmap.end())
						{

							fecesbacmap[species] = Cellfield[i][j][0].volume;
						}
						else
						{
							fecesbacmap[species] += Cellfield[i][j][0].volume;
						}

						Cell Celltemp;
						Cellfield[i][j][k] = Celltemp;
					}
					else
					{
						Cell Celltemp = Cellfield[i][j][k];
						Cellfield[i][j][k] = Cellfield[i + idir][j + jdir][kdir];
						Cellfield[i][j][k].ipos = i;
						Cellfield[i][j][k].jpos = j;
						Cellfield[i + idir][j + jdir][kdir] = Celltemp;
						Cellfield[i + idir][j + jdir][kdir].ipos = i + idir;
						Cellfield[i + idir][j + jdir][kdir].jpos = j + jdir;
					}
					swaps += 1;
					if (swaps == par.maxswaps)
					{
						break;
					}
				}
			}
		}
	}
	for (int i = 0; i < nrow + 2; i++)
	{
		for (int j = 0; j < ncol + 2; j++)
		{
			delete[] hasmoved[i][j];
		}
		delete[] hasmoved[i];
	}
	delete[] hasmoved;
	delete[] randomgrid;
	delete[] randomcell;
	//std::cout<<"Number of swaps performed: " <<swaps <<std::endl;
}

double Grid::Avggrowth() const
{
	/* This method calculates the average growth rate of all bacteria. */
	double growthratetot = 0;
	double numbercells = 0;
	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= ncol; j++)
		{
			for (int k = 0; k < par.maxcells; k++)
			{
				if (Cellfield[i][j][k].exist == true)
				{
					numbercells++;
					growthratetot += Cellfield[i][j][k].growthrate;
				}
			}
		}
	}

	if (numbercells > 0)
		return growthratetot / numbercells;
	else
	{
		cout << "Number of cells equals zero\n";
		return 0;
	}
}

int Grid::Numbercells() const
{
	/* This method calculates the number of bacteria alive in the grid. */
	int numbercells = 0;
	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= ncol; j++)
		{
			for (int k = 0; k < par.maxcells; k++)
			{
				if (Cellfield[i][j][k].exist == true)
				{
					numbercells++;
				}
			}
		}
	}
	return numbercells;
}

void Grid::Record_Host_Uptake(int metabolite, double uptake)
{
	if (uptakemetmap.find(metabolite) == uptakemetmap.end())
	{
		uptakemetmap[metabolite] = -1 * uptake;
	}
	else
	{
		uptakemetmap[metabolite] += -1 * uptake;
	}
}

void Grid::Record_Fecal_Diffusion(int metabolite, double uptake)
{
	for (auto const &x : exchange_metab_dict)
	{
		if (x.second == metabolite)
		{
			if (fecesmetmap.find(x.first) == fecesmetmap.end())
			{
				fecesmetmap[x.first] = uptake;
			}
			else
			{
				fecesmetmap[x.first] += uptake;
			}
		}
	}
}

void Grid::Diffuse_Metabolites(double diffconst)
{
	/* This method calculates the diffusion of metabolites through the grid.
	 * At the borders it also determines the uptake by the host.
	 */
	double ***fieldnew;
	fieldnew = new double **[nrow + 2];
	for (int i = 0; i < nrow + 2; ++i)
	{
		fieldnew[i] = new double *[ncol + 2];
		for (int j = 0; j < ncol + 2; ++j)
		{
			fieldnew[i][j] = new double[numexchangemet];
			for (int k = 0; k < numexchangemet; k++)
			{
				fieldnew[i][j][k] = 0;
			}
		}
	}
	Update_Sumconc();
	double defaultdiff = diffconst;
	double HUrate;
	for (int k = 0; k < numexchangemet; k++)
	{
		if (par.shuffle)
		{
			break;
		}
		HUrate = 0;
		diffconst = defaultdiff;

		for (int i = 1; i <= nrow; i++) // These should be cleared
		{
			metfield[i][0][k] = 0;
			metfield[i][ncol + 1][k] = 0;
		}
		for (int j = 1; j <= ncol; j++) // The ends should also be cleared
		{
			metfield[0][j][k] = 0;
			metfield[nrow + 1][j][k] = 0;
		}
		if (sumconc[k] > par.DIFF_CUTOFF) // Do not diffuse very very rare substances
		{
			for (int i = 0; i <= nrow + 1; i++) // This is the short side
			{
				if (par.shuffle)
				{
					continue; // Don't let metabolites flow out if we're shuffling
				}
				for (int j = 1; j <= ncol; j++) // This is the long side
				{
					fieldnew[i][j][k] = metfield[i][j][k]; // Make a temporary field

					for (int m = 1; m <= 4; m++) // Look seperately for each of the four von neumann neighbours
					{
						if (i == 0 || i == nrow + 1 || j == 0 || j == ncol + 1) // This is the area outside of our actual gut - no diffusion should occur
						{
							continue;
						}
						else if (i == 1 && m == 4)
						{
							fieldnew[i][j][k] += HUrate * diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m)));
							Record_Host_Uptake(k, HUrate * diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m))));
						}
						else if (i == nrow && m == 2)
						{
							fieldnew[i][j][k] += HUrate * diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m)));
							Record_Host_Uptake(k, HUrate * diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m))));
						}
						else if (j == 1 && (m == 1 || m == 2 || m == 4))
						{
							fieldnew[i][j][k] += diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m)));
						}
						else if (j > 1)
						{
							// Here j might be ncol
							if (j == ncol && m == 1)
							{
								Record_Fecal_Diffusion(k, diffconst * (metfield[i][j][k]));
							}
							fieldnew[i][j][k] += diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m)));
							// std::cout<<diffconst * (-metfield[i][j][k] + (*MetFieldNeighbour(i, j, k, m)));
						}
					}
				}
			}
			for (int i = 0; i <= nrow + 1; i++)
			{
				for (int j = 1; j <= ncol; j++)
				{
					metfield[i][j][k] = fieldnew[i][j][k];
				}
			}
		}
	}
	for (int i = 0; i < nrow + 2; i++)
	{
		for (int j = 0; j < ncol + 2; j++)
		{
			delete[] fieldnew[i][j];
		}
		delete[] fieldnew[i];
	}
	delete[] fieldnew;

	int *randrow = new int[nrow];
	int *randcol = new int[ncol];
	for (int n = 0; n < nrow; n++)
	{
		randrow[n] = n;
	}
	for (int n = 0; n < ncol; n++)
	{
		randcol[n] = n;
	}
}
void Grid::DriftCells(QVector<BactParam> bps)
{
	/* This method regulates the drift of the metabolites flowing through the system and
	 * if applicable also the drift of the bacteria.
	 */
	bool celldrift = false;
	if (Uniform() < par.celldrift) // The parameter celldrift determines the whether the bacteria also drift.
		celldrift = true;

	for (int i = 1; i <= nrow; i++)
	{
		if (par.celldrift < 2 || (par.celldrift == 2 && abs(i - (nrow / 2)) < RandNum(nrow / 2)) || (par.celldrift == 3 && Uniform() < (pow(i, 2) / pow(nrow, 2)))) // i> RandNum(nrow+3) && Uniform()>0.33))//slower near sides if celldrift ==2
		{
			for (int j = ncol; j >= 1; j--)
			{

				if ((i > 1 && i < nrow) || (i > 1 && par.celldrift == 3)) // The uppermost and lowermost layer do not drift if celldrift <3, presumably intended to simulate crypts?
				{
					if (celldrift)
					{
						for (int k = 0; k < par.maxcells; k++)
						{
							if (j == ncol && Uniform() < par.cellinflux)
							{
								// par.cellinflux is always 0, the code below should not execute
								if (Cellfield[i][j][k].exist)
								{
									ancestorcount++;
									BactParam bp = bps.at(RandNum(bps.size()) - 1);
									Cellfield[i][1][k] = Cell(bp);
									Cellfield[i][1][k].Init(i, 1);
								}
								else
								{
									Cellfield[i][1][k] = Cell();
								}
								Cellfield[i][j][k] = Cellfield[i][j - 1][k];
								Cellfield[i][j][k].ipos = i;
								Cellfield[i][j][k].jpos = j;
							}
							else
							{
								if (j == (ncol))
								{
									string species = Cellfield[i][j][0].species;
									if (fecesbacmap.find(species) == fecesbacmap.end())
									{
										fecesbacmap[species] = Cellfield[i][j][0].volume;
									}
									else
									{
										fecesbacmap[species] += Cellfield[i][j][0].volume;
									}
								}

								Cellfield[i][j][k] = Cellfield[i][j - 1][k];
								Cellfield[i][j][k].ipos = i;
								Cellfield[i][j][k].jpos = j;
								if (j == 1)
								{
									Cell Celltemp;
									Cellfield[i][j][k] = Celltemp; // No living cells can 'drift' into the gut
								}
							}
						}
					}
				}

				if (j == ncol)
				{
					for (auto const &x : exchange_metab_dict)
					{
						if (x.first != "h2o" && x.first != "o2")
						{
							if (x.second > 0.000001)
							{
								if (fecesmetmap.find(x.first) == fecesmetmap.end())
								{
									fecesmetmap[x.first] = metfield[i][j][x.second];
								}
								else
								{
									fecesmetmap[x.first] += metfield[i][j][x.second];
								}
							}
						}
					}
				}

				for (int k = 0; k < numexchangemet; k++)
				{
					if (k != exchange_metab_dict.at("o2") && k != exchange_metab_dict.at("h2o"))
					{

						metfield[i][j][k] = metfield[i][j - 1][k];
					}
				}
			}
		}
		else
		{
			for (int j = ncol; j >= 1; j--)
			{
				if (j == ncol)
				{
					for (auto const &x : exchange_metab_dict)
					{
						if (x.first != "h2o" && x.first != "o2")
						{
							if (x.second > 0.000001)
							{
								if (fecesmetmap.find(x.first) == fecesmetmap.end())
								{
									fecesmetmap[x.first] = metfield[i][j][x.second];
								}
								else
								{
									fecesmetmap[x.first] += metfield[i][j][x.second];
								}
							}
						}
					}
				}

				for (int k = 0; k < numexchangemet; k++)
				{
					if ((k != exchange_metab_dict.at("o2")) && (k != exchange_metab_dict.at("h2o")))
					{											   // Oxygen gets flushed out much more slowly than actual feces
						metfield[i][j][k] = metfield[i][j - 1][k]; // This moves all metabolites distally by one step
					}
				}
			}
		}
	}
}

void Grid::CellShuffle()
{
	/* This method randomly shuffles the bacteria with a position anywhere in the grid.
	 * If done at every iteration produces a well-mixed situation.
	 * Only used when shuffling, or once with bacplacetype == 3 to distribute bacteria nicely across the field
	 */
	int xa, xb, ya, yb, layera, layerb, pos = 0, target;

	Cell Celltemp;
	long size = ncol * nrow * par.maxcells;
	for (xa = 1; xa <= nrow; xa++)
	{
		for (ya = 1; ya <= ncol; ya++)
		{
			for (layera = 0; layera < par.maxcells; layera++)
			{
				target = pos + RandNum(size) - 1;
				layerb = target % par.maxcells;
				xb = (target) / ncol + 1;
				yb = (target) % ncol + 1; // Calculate coordinates of target
				Celltemp = Cellfield[xa][ya][layera];
				Cellfield[xa][ya][layera] = Cellfield[xb][yb][layerb];
				Cellfield[xb][yb][layerb] = Celltemp;
				Cellfield[xa][ya][layera].ipos = xa;
				Cellfield[xa][ya][layera].jpos = ya;
				Cellfield[xb][yb][layerb].ipos = xb;
				Cellfield[xb][yb][layerb].jpos = yb;
				pos++;
				size--;
			}
		}
	}
}

void Grid::Mix_Metabolites()
{
	/* This method sets the concentrations of all metabolites equal to the average.
	 * If done at every iteration produces a well-mixed situation.
	 * Normally not used!
	 */
	Update_Sumconc();
	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= ncol; j++)
		{
			for (int k = 0; k < numexchangemet; k++)
			{
				metfield[i][j][k] = sumconc[k];
			}
		}
	}
}

void Grid::Set_Unlimited_Metabolite()
{

	// Set some very general compounds to very high levels, these should truly be unlimited
	std::set<string> unlimmetabs = {"h2o"};
	for (std::string m : unlimmetabs)
	{
		if (exchange_metab_dict.count(m) > 0)
		{
			for (int i = 1; i <= nrow; i++)
			{
				for (int j = 1; j <= ncol; j++)
				{
					double localvalue = metfield[i][j][exchange_metab_dict.at(m)];
					if (localvalue < 2777)
					{
						metfield[i][j][exchange_metab_dict.at(m)] = 2777;
						if (inputmetmap.find(m) == inputmetmap.end())
						{
							inputmetmap[m] = 2777 - localvalue;
						}
						else
						{
							inputmetmap[m] += 2777 - localvalue;
						}
					}
				}
			}
		}
	}
}

void Grid::ReadEnv(int i, int j)
{
	/* This method sets the upper bounds of the uptake reactions of the bacteria according to
	 * the locally available concentrations of metabolites.
	 */
	double totcellvol = 0;
	std::set<string> butyrogenics = {"Cbut","Rinu","Ehal"};
	std::set<string> bifids = {"Bbre", "Binf", "Blong"};
	std::set<string> enteros = {"Ecol"};
	std::string thisspecies = Cellfield[i][j][0].species;
	for (int k2 = 0; k2 < par.maxcells; k2++)
	{
		totcellvol += Cellfield[i][j][k2].volume;
	}
	for (int k2 = 0; k2 < par.maxcells; k2++)
	{
		for (auto const &x : exchange_metab_dict)
		{
			int k = Cellfield[i][j][k2].exchange_dict[x.first].reaction_in;
			if (k != 0)
			{
				// Various knockout options, usually disabled
				if (par.knockoutlactate && x.first == "lac_L")
				{
					Cellfield[i][j][k2].ub[k - 1] = 0;
					continue;
				}
				if (par.knockoutnonbiflactose && x.first == "lcts")
				{
					if (bifids.find(thisspecies) == bifids.end())
					{
						Cellfield[i][j][k2].ub[k - 1] = 0;
						continue;
					}
				}
				if (par.knockoutnonbiflactoselate && x.first == "lcts")
				{
					if (simultime > 5040)
					{
						if (bifids.find(thisspecies) == bifids.end())
						{
							Cellfield[i][j][k2].ub[k - 1] = 0;
							continue;
						}
					}
				}
				if (par.knockoutbutyrolactoselate && x.first == "lcts")
				{
					if (simultime > 5040)
					{
						if (butyrogenics.find(thisspecies) != butyrogenics.end()) 
						{
							Cellfield[i][j][k2].ub[k - 1] = 0;
							continue;
						}
					}
				}
				if (par.knockoutenterooxygen && x.first == "o2")
				{
					if (enteros.find(thisspecies) != enteros.end())
					{
						Cellfield[i][j][k2].ub[k - 1] = 0;
						continue;
					}
				}
				if(par.knockout12ppd && x.first == "12ppd_S"){
					Cellfield[i][j][k2].ub[k - 1] = 0;
					continue;
				}

				double maximumuptake = min(par.timescale * metfield[i][j][x.second] / (par.cellfrac * totcellvol), par.max_metup);
				Cellfield[i][j][k2].ub[k - 1] = maximumuptake;
				if (Cellfield[i][j][k2].ub[k - 1] < par.DIFF_CUTOFF)
				{
					Cellfield[i][j][k2].ub[k - 1] = 0;
				}
			}
		}
	}
}

void Grid::HostInflux()
{
	/* This method adds oxygen at the positions at the borders.
	 * How much oxygen is determined by the parameter o2.
	 */
	// Upper and lower border

	if (par.oxin > 0)
	{
		int totalsquares;
		for (int j = 1; j <= ncol; j++)
		{
			// Both long sides are gut wall
			totalsquares = ncol * 2;
			if(par.diminishoxygenrelease){
			double timemod = max(0,1-(simultime/5040));
			metfield[1][j][exchange_metab_dict.at("o2")] += (par.oxin*timemod) / totalsquares;
			metfield[nrow][j][exchange_metab_dict.at("o2")] += (par.oxin*timemod) / totalsquares;
			}
			else{
			metfield[1][j][exchange_metab_dict.at("o2")] += par.oxin / totalsquares;
			metfield[nrow][j][exchange_metab_dict.at("o2")] += par.oxin / totalsquares;

			}
		}
		// left border, i.e. small intestine output. Could also contain oxygen
		for (int i = 1; i <= nrow; i++)
		{
			// metfield[i][1][exchange_metab_dict.at("o2")] += par.o2/totalsquares;
		}
	}
};
void Grid::MixedSugarInflux()
{
	/* This method adds sugar to the left end of the grid
	 */

	double intakemodifier; // Can be used to scale the milk quantity consumed with age
	intakemodifier = par.baseintake + (par.intakeincrease * (simultime - 1));
	if (intakemodifier < 0)
	{
		intakemodifier = 0;
	}
	std::cout << "Current time is " << simultime << std::endl;

	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= 6; j++) // Release food only in the first 6 squares
		{
			if (par.foodinglucose > 0)
			{
				metfield[i][j][exchange_metab_dict.at("glc_D")] += par.foodinglucose / (6 * nrow) * intakemodifier;
				if (inputmetmap.find("glc_D") == inputmetmap.end())
				{
					inputmetmap["glc_D"] = par.foodinglucose / (6 * nrow) * intakemodifier;
				}
				else
				{
					inputmetmap["glc_D"] += par.foodinglucose / (6 * nrow) * intakemodifier;
				}
			}
			if (par.foodingalactose > 0)
			{
				metfield[i][j][exchange_metab_dict.at("gal")] += par.foodingalactose / (6 * nrow) * intakemodifier;
				if (inputmetmap.find("gal") == inputmetmap.end())
				{
					inputmetmap["gal"] = par.foodingalactose / (6 * nrow) * intakemodifier;
				}
				else
				{
					inputmetmap["gal"] += par.foodingalactose / (6 * nrow) * intakemodifier;
				}
			}

			if (par.foodinlactose > 0)
			{
				metfield[i][j][exchange_metab_dict.at("lcts")] += par.foodinlactose / (6 * nrow) * intakemodifier;
				if (inputmetmap.find("lcts") == inputmetmap.end())
				{
					inputmetmap["lcts"] = par.foodinlactose / (6 * nrow) * intakemodifier;
				}
				else
				{
					inputmetmap["lcts"] += par.foodinlactose / (6 * nrow) * intakemodifier;
				}
			}

			if (par.foodinlactate > 0)
			{
				metfield[i][j][exchange_metab_dict.at("lac_L")] += par.foodinlactate / (6 * nrow) * intakemodifier;
				if (inputmetmap.find("lac_L") == inputmetmap.end())
				{
					inputmetmap["lac_L"] = par.foodinlactate / (6 * nrow) * intakemodifier;
				}
				else
				{
					inputmetmap["lac_L"] += par.foodinlactate / (6 * nrow) * intakemodifier;
				}
			}
			if (par.foodinfl > 0)
			{
				metfield[i][j][exchange_metab_dict.at("2FuLa")] += par.foodinfl / (6 * nrow) * intakemodifier;
				if (inputmetmap.find("2FuLa") == inputmetmap.end())
				{inputmetmap["2FuLa"] =  par.foodinfl / (6 * nrow) * intakemodifier;}
				else
				{inputmetmap["2FuLa"] +=  par.foodinfl / (6 * nrow) * intakemodifier;}
			}

			if (par.foodinGOS> 0)
			{
				double gospertile = (double)((par.foodinGOS) / (double)(6 * nrow)) * intakemodifier; //no *100 factor here
				metfield[i][j][exchange_metab_dict.at("GOS5")] += gospertile / 12.267;
				metfield[i][j][exchange_metab_dict.at("GOS4")] += gospertile / 3.608;
				metfield[i][j][exchange_metab_dict.at("GOS3")] += gospertile / 1.559;
				
				if (inputmetmap.find("GOS3") == inputmetmap.end())
				{inputmetmap["GOS3"] =  gospertile / 1.559;}
				else
				{inputmetmap["GOS3"] +=  gospertile / 1.559;}
				if (inputmetmap.find("GOS4") == inputmetmap.end())
				{inputmetmap["GOS4"] =  gospertile / 3.608;}
				else
				{inputmetmap["GOS4"] +=  gospertile / 3.608;}
				if (inputmetmap.find("GOS5") == inputmetmap.end())
				{inputmetmap["GOS5"] =  gospertile / 12.267;}
				else
				{inputmetmap["GOS5"] +=  gospertile / 12.267;}
			}
		}
	}
}

void Grid::Update_Concentrations(int i, int j)
{
	/* This method updates the local concentrations according to what the bacteria at that position take up and excrete.
	 */

	for (int k = 0; k < par.maxcells; k++)
	{
		for (auto const &x : exchange_metab_dict)
		{
			// Input concentrations not yet realistic and measured in: mmol/g dcw.
			double foodup = 0;
			int r_in = Cellfield[i][j][k].exchange_dict[x.first].reaction_in;
			int r_out = Cellfield[i][j][k].exchange_dict[x.first].reaction_out;
			if (r_in != 0 && r_out != 0)
			{
				foodup = (Cellfield[i][j][k].flux[r_in - 1] - Cellfield[i][j][k].flux[r_out - 1]) * (Cellfield[i][j][k].volume * par.cellfrac / par.timescale);
			}
			else if (r_in != 0)
			{
				foodup = Cellfield[i][j][k].flux[r_in - 1] * (Cellfield[i][j][k].volume * par.cellfrac / par.timescale);
			}
			else if (r_out != 0)
			{
				foodup = (-1.0 * Cellfield[i][j][k].flux[r_out - 1]) * (Cellfield[i][j][k].volume * par.cellfrac / par.timescale);
			}
			metfield[i][j][x.second] -= foodup;
			if (metfield[i][j][x.second] < 0) // due to floating point errors in the linear solver
			{
				if (metfield[i][j][x.second] < -0.00000000001)
				{
					std::cout << "Negative metabolite recorded" << std::endl;
					std::cout << "It's " << x.first << " at " << metfield[i][j][x.second] << std::endl;
					ofstream errorfile;
					errorfile.open(simulationdir + "/errorlog.txt", fstream::app);
					errorfile << "Negativemetaberror"
							  << "," << Cellfield[i][j][k].species << "," << x.first << "," << metfield[i][j][x.second] << "," << simultime << "," << i << "," << j << std::endl;
					errorfile.close();
				}
				metfield[i][j][x.second] = 0;
			}
		}
	}
}

void Grid::Step_Metabolism(QVector<BactParam> bps)
{
	/* Let all bacteria in the entire grid metabolize and update the concentrations of metabolites accordingly.
	 */
	int *randomgrid = new int[nrow * ncol];
	int *randomcell = new int[par.maxcells];
	for (int i = 0; i < nrow * ncol; i++)
	{
		randomgrid[i] = i;
	}
	for (int i = 0; i < par.maxcells; i++)
	{
		randomcell[i] = i;
	}
	Shufflevector(randomgrid, nrow * ncol);

	for (int pos = 0; pos <= nrow * ncol - 1; pos++)
	{
		// field-loop
		int i = randomgrid[pos] / (ncol) + 1;
		int j = randomgrid[pos] % (ncol) + 1;

		// FBA metabolism
		Shufflevector(randomcell, par.maxcells); // access the layers in random order
		for (int k = 0; k < par.maxcells; k++)
		{
			int kcell = randomcell[k];
			if (Cellfield[i][j][kcell].exist == true) // if cell exist in this layer
			{
				Cellfield[i][j][kcell].active = false;
				Cellfield[i][j][kcell].time_alive++;
				ReadEnv(i, j); // read external environment for fba
				if (Cellfield[i][j][kcell].volume < par.divisionsize * 2)
				{
					Cellfield[i][j][kcell].active = true;
				}
				if (Cellfield[i][j][kcell].active == true)
				{
					Cellfield[i][j][kcell].Do_Metabolism(bps);
				}
				else // if Cell is not active
				{
					Cellfield[i][j][kcell].growthrate = 0;
					for (int k2 = 0; k2 < Cellfield[i][j][kcell].numreactions; k2++)
					{
						Cellfield[i][j][kcell].flux[k2] = 0;
					}
				}
			}
		}

		Update_Concentrations(i, j); // This applies the FBA result (the fluxes) to the environment
		for (int k = 0; k < par.maxcells; k++)
		{
			if (Cellfield[i][j][k].exist == true) // if cell exist in this layer
			{
				Cellfield[i][j][k].volume += Cellfield[i][j][k].volume * Cellfield[i][j][k].growthrate / par.timescale;
			}
		} // loop over layers

	} // end of loop over field

	delete[] randomgrid;
	delete[] randomcell;
}

void Grid::Grid_Die_Reproduce()
{
	/* This method regulates density dependent cell death and division of the bacteria.
	 */
	int *randomgrid = new int[nrow * ncol];
	for (int i = 0; i < nrow * ncol; i++)
	{
		randomgrid[i] = i;
	}

	Shufflevector(randomgrid, nrow * ncol);

	for (int pos = 0; pos <= nrow * ncol - 1; pos++)
	{
		int i = randomgrid[pos] / (ncol) + 1;
		int j = randomgrid[pos] % (ncol) + 1;
		// random Cell death
		for (int k = 0; k < par.maxcells; k++)
		{
			if (Cellfield[i][j][k].exist == true && Uniform() < ((par.pdeath_base / par.timescale) + CountCellNeighbour(i, j) * (par.pdeath / par.timescale)))
			{
				Cellfield[i][j][k].Die();
			}
		}
		while (Cellfield[i][j][0].volume > par.divisionsize && CountCellNeighbour(i, j) < 1.0)
		{ // While volume large, and there is at least one empty spot
			// Because it's a while loop, we'll keep trying until succesfull
			//  This part only works for one bacterium per position at the moment.
			//  If in the future more bacteria per position are used this part needs to be adjusted.
			int a = RandNum(3) - 2;
			int b = RandNum(3) - 2;
			if ((i + a) >= 1 && (i + a) <= nrow && (j + b) >= 1 && (j + b) <= ncol)
			{
				if (Cellfield[i + a][j + b][0].exist == false)
				{
					double v1 = Cellfield[i][j][0].volume * 0.5; // JUNE2020 * (0.3 + (0.4 * Uniform()));
					double v2 = Cellfield[i][j][0].volume - v1;
					Cellfield[i + a][j + b][0] = Cellfield[i][j][0];
					Cellfield[i][j][0].volume = v1;
					Cellfield[i + a][j + b][0].volume = v2;
					Cellfield[i + a][j + b][0].numancestors += 1;
				}
			}
		}
	}

	delete[] randomgrid;
}

void Grid::Place_Bacterium(QVector<BactParam> bps)
{ // Places bacteria, location depending on immigration parameters
	if (par.immigration > 0)
	{
		if (par.immigration > 1 || simultime < 1000) // Stop placing after 1000 steps unless par.immigration is high
		{
			if (par.immigration > 2)
			{
				if (par.immigration < 4 || simultime < 6720)
				{
					for (int i = 0; i <= nrow; i++)
					{
						for (int j = 0; j < ncol; j++)
						{
							if (!(Cellfield[i][j][0].indexnr > 0)) // Can't place if field is in use
							{
								if (Uniform() < par.placechance)
								{
									BactParam bp;

									if (par.bacplacetype == 3)
									{
										// This could be more efficient

										std::string thisspecies = bacspecies[RandNum(bacspecies.size()) - 1];
										for (int b = 0; b < bps.length(); b++)
										{
											if (bps.at(b).speciesName == thisspecies)
											{
												bp = bps.at(b);
											}
										}
									}
									else
									{

										bp = bps.at(RandNum(bps.size()) - 1); // Place random bacterium
									}
									// BactParam bp = bps.at(bactindex); //Place specific bacterium
									std::cout << "Placed " << bp.speciesName << " at X=" << i << " Y=" << j << std::endl;
									Cellfield[i][j][0] = Cell(bp);
									Cellfield[i][j][0].Init(i, j);
								}
							}
						}
					}
				}
			}
			else
			{
				for (int i = 0; i <= nrow; i++)
				{
					// std::cout << Cellfield[i][1][0].indexnr<< std::endl;
					if (Uniform() < par.placechance)
					{
						if (!(Cellfield[i][1][0].indexnr > 0))
						{
							BactParam bp = bps.at(RandNum(bps.size()) - 1); // Place random bacterium
							// BactParam bp = bps.at(bactindex); //Place specific bacterium
							// std::cout << "Placed " << bp.speciesName << " at X="<< i << " Y="<<1<< std::endl;
							Cellfield[i][1][0] = Cell(bp);
							Cellfield[i][1][0].Init(i, 1);
						}
					}
				}
			}
		}
	}
}

void Grid::UsedReactions()
{
	/* Print which reactions of their metabolic network the bacteria in the grid actually use.
	 * Used for testing and such.
	 */
	for (int i = 1; i <= nrow; i++)
	{
		for (int j = 1; j <= ncol; j++)
		{
			std::cout << Cellfield[i][j][0].species << " " << i << " " << j << " " << std::endl;
			for (auto const &x : Cellfield[i][j][0].react_dict)
			{
				if (Cellfield[i][j][0].used[x.second - 1] == 1)
				{
					std::cout << x.first << std::endl;
				}
			}
		}
	}
}
