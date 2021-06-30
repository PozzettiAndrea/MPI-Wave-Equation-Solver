#include "Funcs.h"

using namespace std;

//PROBLEM VARIABLES//
double **grid, **new_grid, **old_grid;         		//Three grids, one for each time level
int imax = 100, jmax = 100;                    		//Rows and columns in the global domain
double t_max = 30.0;						   		//Maximum time to reach with iterations
double t, t_out = 0.0, dt_out = 0.25, dt;
double y_max = 10.0, x_max = 10.0, dx, dy;	   		//Actual dimensions of the domain, spatial steps
double c = 1;                                  		//Wave speed

int id, p;									   		//Core id and total number of cores

int drows, dcols;                  			   		//Global domain rows and columns
int sdrows, sdcols;                 		   		//Subdomain rows and columns
int global_i, global_j;         			   		//Coordinates of the global origin of the subdomain
int sd_i, sd_j;                      		   		//Subdomain row and column index

int nb_trbl[4];   							   		//Neighbour ids: top, right, bottom and left
bool boundary_trbl[4] = {false,false,false,false};	//Boundary situation: is the top, right, bottom or left side touching a (global) boundary?
bool communicating[4] = {true, true, true, true};	//Communication situation: should we start a comm for our top, right, bottom or left sides?
std::vector<MPI_Request> request;          	   		//MPI Request vector
MPI_Datatype l_send, l_recv, r_send, r_recv, t_send, t_recv, b_send, b_recv;  //Datatypes to send and receive our local boundaries

int btype;                               	   		//Setting Periodic, Dirichlet or Neumann BCs

int nsources;                				   		//Number of sources
std::vector<std::vector<double>> sources;	   		//Vector to store global sources
std::vector<std::vector<double>> local_sources;		//Vector to store local sources

bool printing = true;                          		//Switch for the HPC

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

	double start = MPI_Wtime();
	
	//----READING OUR INPUT----//


	//Reading imax, jmax, boundary type, sources
	try
	{	
		imax = std::atoi(argv[1]);
		jmax = std::atoi(argv[2]);
		btype = std::atoi(argv[3]);
		nsources = std::atoi(argv[4]);
		
		if (imax < 1 || jmax <1  )
		{
			std::cout<<"imax and jmax cannot be negative values!\n";
			throw "imaxjmax";
		}

		if (btype != 0 && btype != 1 && btype != 2)
		{
			std::cout<<"Sorry, btype should be an integer between 0 and 2 (0 for Parallel, 1 for Dirichlet, 2 for Neumann)\n";
			throw "btype";
		}

		if (nsources < 1)
		{
			std::cout<<"Number of sources must be at least 1\n";
			throw "nsources";
		}
			
	} catch(...){
		std::cout << "Please repeat your input. Format is: imax jmax btype nsources\n";
		return 0;
	}



	//Reading the sources
	for (int i = 0; i<nsources;i++)
	{
		std::vector<double> source;
		try
		{
			source.push_back(stod(argv[4+5*i+1]));
			source.push_back(stod(argv[4+5*i+2]));
			source.push_back(stod(argv[4+5*i+3]));
			source.push_back(stod(argv[4+5*i+4]));
			source.push_back(stod(argv[4+5*i+5]));
		} catch (...)
		{
			std::cout<< "Error reading source list\n";
			return 0;
		}

		if (source.size() == 5)
		{
			sources.push_back(source);
		} else {
			std::cout<< "Error: Please provide sources in correct format: i j period intensity time\n";
			return 0;
		}
	}
	
	//----READING OUR INPUT----//

	

	//----SETTING UP THE PROBLEM VARIABLES----//

	//Finding the dimensions of our subdivision grid
	find_dimensions(p, drows, dcols);

	//cout << "Rows: " << drows << " Columns: " << dcols << "\n";

	//Finding our subdomain's i and j indices within the subdivision grid
	sd_i = id / dcols;
    sd_j = id % dcols;

	neighbours();

	//Finding our subdomain rows, columns and global origin coordinates
	sdrows = ((sd_i + 1) * imax / drows) - (sd_i * imax / drows);
    sdcols = ((sd_j + 1) * jmax / dcols) - (sd_j *jmax / dcols);
	global_i = (sd_i * imax / drows);
	global_j = (sd_j * jmax / dcols);
	
	//We allocate the space for our pointer pointer grids
	fill_ppointer(old_grid);
	fill_ppointer(grid);
	fill_ppointer(new_grid);

	//We find out the spatial step sizes
	dx = x_max / ((double)sdrows - 1);
	dy = y_max / ((double)sdcols - 1);

	//Time set to 0
	t = 0.0;

	//The timestep
	dt = 0.1 * min(dx, dy) / c;

	//Printing our first grid
	int out_cnt = 0, it = 0;
	grid_to_file(out_cnt);
	out_cnt++;
	t_out += dt_out;

	//Are we touching any global boundaries?
    BoundaryCheck();
	
	//Finding out which of the global sources are local
	set_sources();
	
	//----SETTING UP THE PROBLEM VARIABLES----//



	//----MAIN LOOP----//
	while (t < t_max)
	{
		do_iteration();

		if (printing && t_out <= t)
		{
			// if (id ==0)
			// {
			// 	std::cout << "output: " << out_cnt << "\tt: " << t << "\titeration: " << it << endl;
			// }

			grid_to_file(out_cnt);
			out_cnt++;
			t_out += dt_out;
		}

		it++;
	}
	//----MAIN LOOP----//


	double end = MPI_Wtime();

	//Couting if id == 0
	if (id == 0)
	{
		std::cout<< "Took: " <<  end-start << "s for a " << imax << " x " << jmax << "grid with " << p << " cores\n";
	}

	MPI_Finalize();

	//----FREEING MEMORY----//
    del_ppointer(old_grid);
	del_ppointer(grid);
	del_ppointer(new_grid);
	//----FREEING MEMORY----//
}

