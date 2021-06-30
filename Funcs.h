#include <iostream>
#include <vector>
#include <algorithm> 
#include <sstream>
#include <fstream>
#include <cmath>
#include <mpi.h>


extern double **grid, **new_grid, **old_grid;           //Three grids, one for each time level
extern int imax, jmax;                                  //Rows and columns in the global domain
extern double t_max;                                    //Maximum time to reach with iterations
extern double t, t_out, dt_out, dt;  
extern double y_max, x_max, dx, dy;                     //Actual dimensions of the domain, spatial steps
extern double c;                                        //Wave speed

extern int id, p;                                       //Core id and total number of cores

extern int drows, dcols;                                //Global domain rows and columns
extern int sdrows, sdcols;                              //Subdomain rows and columns
extern int global_i, global_j;                          //Coordinates of the global origin of the subdomain
extern int sd_i, sd_j;                                  //Subdomain row and column index

extern int nb_trbl[4];                                  //Neighbour ids: top, right, bottom and left
extern bool boundary_trbl[4];                           //Boundary situation: is the top, right, bottom or left side touching a (global) boundary?
extern bool communicating[4];                           //Communication situation: should we start a comm for our top, right, bottom or left sides?
extern std::vector<MPI_Request> request;                //MPI Request vector
extern MPI_Datatype l_send, l_recv, r_send, r_recv, t_send, t_recv, b_send, b_recv; //Datatypes to send and receive our local boundaries

extern int btype;                                       //Setting Periodic, Dirichlet or Neumann BCs

extern int nsources;                                    //Number of sources
extern std::vector<std::vector<double>> sources;        //Vector to store global sources
extern std::vector<std::vector<double>> local_sources;  //Vector to store local sources

//Finding the shape of our subdivision grid
void find_dimensions(int p, int &drows, int &dcols);

//Filling a pointer with a grid
void fill_ppointer(double **&ptr);

//Deallocating the memory from the grid
void del_ppointer(double **&ptr);

//Printing to file
void grid_to_file(int out);

//Check if we are close to any boundary
void BoundaryCheck();

//Find the ids of the neighbours
void neighbours();

//Finding core id from global subdivision indices
int toid(int sdi, int sdj);

//Iterate once
void do_iteration();

//Apply boundary conditions
void applybcs();

//Dirichlet boundary conditions
void Dirichlet();

//Neumann boundary conditions
void Neumann();

//Paralle boundary conditions
void Parallel();

//Create MPI types for boundaries
void create_types();

// Do the round of communications after iterating.
void communicate();

//Set the sources: iterate through all sources to find the ones in the local domain
void set_sources();

//Implement the source
void source(int source_i, int source_j, double period, double intensity, double time);