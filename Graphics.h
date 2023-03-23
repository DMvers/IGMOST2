#ifndef GRAPHICS_H
#define GRAPHICS_H
#include "Grid.h"

/* The Graphics class produces images that visualise the grid (the bacteria and
 * the concentrations of different metabolites), these images can be turned into a video.
 */

class Graphics
{
private:
    int ***fields_display; // array with grids for each visualization
    int countmovie;        // the number of the current frame/iteration
    int numfields;         // number of visualizations (one for the bacteria and one for each metabolite)
    int nrow;              // y-dimension of the grid
    int ncol;              // x-dimension of the grid
    std::string moviedir;  // the name of the movie directory
    Cell *Celllayer;       // pointer to cell that is currently visualized

public:
    Graphics(const Grid &g);                // This method constructs the Graphics class on the basis of an instance of the Grid.
    Graphics(const Graphics &g);            // This method constructs the Graphics class on the basis of a previous instance of Graphics.
    ~Graphics();                            // The destructor of the Graphics class.
    Graphics &operator=(const Graphics &g); // Making a reference of the Graphics class.
    void Init_Movie();                      // Initialisation of the Graphics class. The movie directory is created.
    void Create_Planes(const Grid &g);      // This method creates the different visualisations of the grid for the movie.
    void Write_Movie(const Grid &g);        // This method renders the actual PNGs for the movie.
};

#endif
