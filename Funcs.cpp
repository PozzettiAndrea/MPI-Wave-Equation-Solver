#include "Funcs.h"

//Finding the shape of our subdivision grid
void find_dimensions(int p, int& drows, int& dcols)
{
    //Initial "dummy" values
	int min_gap = imax;
    int gap = imax;

    //Top integer to check
	int top = sqrt(p) + 1;

    //Check if any of the integers up to "top" is a divisor of p
	for (int i = 1; i <= top; i++)
	{
		if (p % i == 0)
		{
            //We are seeking to minimise the difference in length (as measured in cells) of subdomain boundaries:
            //If a domain is 600x100 and we have 6 cores, the optimal solution is having a long row of 6 cores

			//Checking vertical orientation (2x3)
            gap = abs((jmax/(p / i)) - imax/i);

			if (gap < min_gap)
			{
				min_gap = gap;
				drows = i;
                dcols = p / i;
			}

            //Checking horizontal orientation (3x2)
            gap = abs((jmax/(i)) - imax/(p / i));

			if (gap < min_gap)
			{
				min_gap = gap;
				drows = p / i;
                dcols = i;
			}
		}
	}
}

//Filling a pointer with a grid
void fill_ppointer(double **&ptr)
{
    // std::cout << sdrows << " " << sdcols <<"\n";

    //Fill the pointer pointer with pointers
	ptr = new double*[sdrows+2];
    
	for (int i = 0; i < sdrows+2; i++)
	{
		ptr[i] = new double[sdcols+2];
        std::fill_n(ptr[i], sdcols+2, 0);
	}
}

//Deallocating the memory from the grid
void del_ppointer(double **&ptr)
{
	for (int i = 0; i < sdrows+2; i++)
	{
        delete[] ptr[i];
    }
	delete[] ptr;
}

//Printing to file
void grid_to_file(int out)
{
	std::stringstream fname;
	std::fstream f1;

    //We write the rows and the columns in each file name
	fname << "./out/id" << "_" << id << "_" << drows << "x" << dcols << "_output_" << out << ".dat";

    f1.open(fname.str().c_str(), std::ios_base::out);
	for (int i = 1; i < sdrows + 1; i++)
	{
		for (int j = 1; j < sdcols + 1; j++)
			f1 << grid[i][j] << "\t";
		f1 << std::endl;
	}
	f1.close();
}

//Check if we are close to any boundary
void BoundaryCheck()
{
    if (sd_i == 0) //Top row
    {
        boundary_trbl[0] = true;
    }
    if (sd_j == dcols - 1) //Right edge
    {
        boundary_trbl[1] = true;
    }
    if (sd_i == drows - 1) //Bottom row 
    {
        boundary_trbl[2] = true;
    }
    if (sd_j == 0) //Left edge
    {
        boundary_trbl[3] = true;
    }

}

//Find the ids of the neighbours
void neighbours()
{
    nb_trbl[0] = toid(sd_i - 1, sd_j);
	nb_trbl[1] = toid(sd_i, sd_j + 1);
	nb_trbl[2] = toid(sd_i + 1, sd_j);
	nb_trbl[3] = toid(sd_i, sd_j - 1);
}

//Finding core id from global subdivision indices
int toid(int sdi, int sdj)
{
    if (btype != 0) //If the boundaries aren't parallel, don't try out-of-bounds indices    
    {
        if(sdi < 0 || sdj < 0 || sdi-1 > drows || sdj-1 > dcols)
        {
            return -1;
        }
        else
        {
            return dcols*sdi + sdj;
        }
    }
    else //if boundaries are parallel, just add drows/dcols and take the mod
    {
        return dcols*((sdi + drows)%drows) + (sdj + dcols)%dcols;
    }
}


//Iterate once
void do_iteration()
{
    //Check local sources
    for (int i = 0; i < local_sources.size(); i++)
    {
        source(local_sources[i][0], local_sources[i][1], local_sources[i][2], local_sources[i][3], local_sources[i][4]);
    }

    //Update the new grid
    //pow() is gone as suggested
	for (int i = 1; i < sdrows + 1; i++)
		for (int j = 1; j < sdcols + 1; j++)
            new_grid[i][j] = (dt * c)*(dt * c) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy)) + 2.0 * grid[i][j] - old_grid[i][j];

    //Apply boundary conditions
    applybcs();

    //Create MPI Datatypes to send and receive
    create_types();

    //Communicate between subdomains
    communicate();

    //Update time
    t += dt;

    //Swap grid indices
    std::swap(old_grid, new_grid);
    std::swap(old_grid, grid);

    
}


//Apply boundary conditions
void applybcs()
{
    if (btype == 0)
    {
        Parallel();
    }

    if (btype == 1)
    {
        Dirichlet();
    } 

    if (btype == 2)
    {
        Neumann();
    }   
    
}


//Dirichlet boundary conditions
void Dirichlet()
{
    if (boundary_trbl[0])
    {
        for (int j = 1; j - 1 < sdcols; j++)
        {
            new_grid[1][j] = 0;
        }

        //Shouldn't try to talk to a wall
        communicating[0] = false;
    }

    if (boundary_trbl[1])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][sdcols] = 0;
        }

        //Shouldn't try to talk to a wall
        communicating[1] = false;
    }

    if (boundary_trbl[2])
    {
        for (int j = 1; j - 1 < sdcols; j++)
        {
            new_grid[sdrows][j] = 0;
        }
        
        //Shouldn't try to talk to a wall
        communicating[2] = false;
    }

    if (boundary_trbl[3])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][1] = 0;
        }

        //Shouldn't try to talk to a wall
        communicating[3] = false;
    }
}

//Neumann boundary conditions
void Neumann()
{
    if (boundary_trbl[0])
    {
        for (int j = 1; j - 1 < sdcols; j ++)
        {
            new_grid[1][j] = new_grid[2][j];
        }

        //Shouldn't try to talk to a wall
        communicating[0] = false;
    }

    if (boundary_trbl[1])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][sdcols - 1] = new_grid[i][sdcols];
        }

        //Shouldn't try to talk to a wall
        communicating[1] = false;
    }

    if (boundary_trbl[2])
    {
        for (int j = 1; j - 1 < sdcols; j++)
        {
            new_grid[sdrows][j] = new_grid[sdrows - 1][j];
        }

        //Shouldn't try to talk to a wall
        communicating[2] = false;
    }

    if (boundary_trbl[3])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][1] = new_grid[i][2];
        }

        //Shouldn't try to talk to a wall
        communicating[3] = false;
    }
}

//Parallel boundary conditions
void Parallel()
{
    if (id == nb_trbl[0])
    {
        for (int j = 1; j - 1 < sdcols; j++)
        {
            new_grid[0][j] = new_grid[sdrows + 1][j];
        }

        //Talking to ourselves already
        communicating[0] = false;
    }

    if (id == nb_trbl[1])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][sdcols + 1] = new_grid[i][1];
        }

        //Talking to ourselves already
        communicating[1] = false;
    }
    
    if (id == nb_trbl[2])
    {
        for (int j = 1; j - 1< sdcols; j++)
        {
            new_grid[sdrows + 1][j] = new_grid[1][j];
        }

        //Talking to ourselves already
        communicating[2] = false;
    }
    
    if (id == nb_trbl[3])
    {
        for (int i = 1; i - 1 < sdrows; i++)
        {
            new_grid[i][0] = new_grid[i][sdcols + 1];
        }

        //Talking to ourselves already
        communicating[3] = false;
    }   
}

//Create MPI types for boundaries
void create_types()
{
    //Create start
    MPI_Aint add_start;
    MPI_Get_address(new_grid, &add_start);
	MPI_Datatype typeval = MPI_DOUBLE;

	MPI_Aint address;
    
    //Creating the send type for the top boundary
	MPI_Get_address(&new_grid[1][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &sdcols, &address, &typeval, &t_send);
	MPI_Type_commit(&t_send);

    //Creating the receive type for the top boundary
	MPI_Get_address(&new_grid[0][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &sdcols, &address, &typeval, &t_recv);
	MPI_Type_commit(&t_recv);

    //Creating the send type for the bottom boundary
    MPI_Get_address(&new_grid[sdrows][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &sdcols, &address, &typeval, &b_send);
	MPI_Type_commit(&b_send);

    //Creating the receive type for the bottom boundary
    MPI_Get_address(&new_grid[sdrows + 1][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &sdcols, &address, &typeval, &b_recv);
	MPI_Type_commit(&b_recv);

    int block_lengths[sdrows];
    MPI_Datatype typelist[sdrows];

    std::fill_n(block_lengths, sdrows, 1);
    std::fill_n(typelist, sdrows, MPI_DOUBLE);

    MPI_Aint sr[sdrows];
    MPI_Aint rr[sdrows];
    MPI_Aint sl[sdrows];
    MPI_Aint rl[sdrows];

    //Send and receive for right and left boundaries
    for (int i = 1; i < sdrows + 1; i++)
    {
		MPI_Aint temp_r;
        MPI_Aint temp_l;
		MPI_Get_address(&new_grid[i][sdcols], &temp_r);
        MPI_Get_address(&new_grid[i][1], &temp_l);
        sr[i-1] = temp_r - add_start;
        rr[i-1] = temp_r + sizeof(double) - add_start; //adding size of double should take us to the next value in the row
        sl[i-1] = temp_l - add_start;
        rl[i-1] = temp_l - 8 - add_start;
    }

	MPI_Type_create_struct(sdrows, &block_lengths[0], &sr[0], &typelist[0], &r_send);
	MPI_Type_create_struct(sdrows, &block_lengths[0], &rr[0], &typelist[0], &r_recv);
    MPI_Type_create_struct(sdrows, &block_lengths[0], &sl[0], &typelist[0], &l_send);
	MPI_Type_create_struct(sdrows, &block_lengths[0], &rl[0], &typelist[0], &l_recv);
    MPI_Type_commit(&r_send);
	MPI_Type_commit(&r_recv);
    MPI_Type_commit(&l_send);
    MPI_Type_commit(&l_recv);

}

// Do the round of communications after iterating.
void communicate()
{
    //Create all the communications, receives before sends :)
    if (communicating[0])
    {
        MPI_Request tempr;
        MPI_Irecv(new_grid, 1, t_recv, nb_trbl[0], 2, MPI_COMM_WORLD, &tempr);
        request.push_back(tempr);

        MPI_Request temps;
        MPI_Isend(new_grid, 1, t_send, nb_trbl[0], 0, MPI_COMM_WORLD, &temps);
        request.push_back(temps);
    }
    if (communicating[1])
    {
        MPI_Request tempr;
        MPI_Irecv(new_grid, 1, r_recv, nb_trbl[1], 3, MPI_COMM_WORLD, &tempr);
        request.push_back(tempr);

        MPI_Request temps;
        MPI_Isend(new_grid, 1, r_send, nb_trbl[1], 1, MPI_COMM_WORLD, &temps);
        request.push_back(temps);
    }
    if (communicating[2])
    {
        MPI_Request tempr;
        MPI_Irecv(new_grid, 1, b_recv, nb_trbl[2], 0, MPI_COMM_WORLD, &tempr);
        request.push_back(tempr);

        MPI_Request temps;
        MPI_Isend(new_grid, 1, b_send, nb_trbl[2], 2, MPI_COMM_WORLD, &temps);
        request.push_back(temps);

    }
    if (communicating[3])
    {
        MPI_Request tempr;
        MPI_Irecv(new_grid, 1, l_recv, nb_trbl[3], 1, MPI_COMM_WORLD, &tempr);
        request.push_back(tempr);

        MPI_Request temps;
        MPI_Isend(new_grid, 1, l_send, nb_trbl[3], 3, MPI_COMM_WORLD, &temps);
        request.push_back(temps);
    }

    //Wait until they have all finished
    MPI_Waitall(request.size(), request.data(), MPI_STATUS_IGNORE);

    //Clear the request array for later use
    request.clear();

    //Freeing the MPI Datatypes for efficiency
    MPI_Type_free(&t_send);
    MPI_Type_free(&t_recv);
    MPI_Type_free(&r_send);
    MPI_Type_free(&r_recv);
    MPI_Type_free(&b_send);
    MPI_Type_free(&b_recv);
    MPI_Type_free(&l_send);
    MPI_Type_free(&l_recv);
}


//Set the sources: iterate through all sources to find the ones in the local domain
void set_sources()
{
    for (int i=nsources; i>0 ;i--)
    {
        if (sources[i-1][0] >= global_i && sources[i-1][0] - global_i < sdrows)
        {
            if (sources[i-1][1] >= global_j && sources[i-1][1] - global_j < sdcols)
            {
                std::vector<double> source;
                source.push_back(1 + sources[i-1][0] - global_i);
                source.push_back(1 + sources[i-1][1] - global_j);
                source.push_back(sources[i-1][2]);
                source.push_back(sources[i-1][3]);
                source.push_back(sources[i-1][4]);
                local_sources.push_back(source);
            }
        }
    }
}

//Implement the source
void source(int source_i, int source_j, double period, double intensity, double time)
{
    if (t <= time)
    {
        grid[source_i][source_j] = intensity * std::cos(6.28*t/period);
        old_grid[source_i][source_j] = intensity * std::cos(6.28*t/period);
    }
}