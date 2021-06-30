## void  find_dimensions(int  p, int  &drows, int  &dcols);
Finding the shape of our subdivision grid. You input p (the number of cores) and the number of rows and the number of columns is written in drows and dcols' address.
It seeks to minimise the difference between each subdomain's number of rows and columns
## void  fill_ppointer(double  **&ptr);
This function allocates space to the grid **ptr passed to it.
I thought it would massively speed up communications from STL but I'm not sure it did.
## void  del_ppointer(double  **&ptr);
This function deallocates the space allocated to the **ptr passed to it.
## void  grid_to_file(int  out);
This function outputs the grid to a .dat file.
The int out is one of the parameters used for choosing the filename
## void  BoundaryCheck();
This function checks if we are close to any global boundary.
If so, it changes a value in bool boundary_trbl[4]
## void  neighbours();
This function checks the ids of the subdomain's neighbours.
It stores them in int nb_trbl[4]
## int  toid(int  sdi, int  sdj);
This function returns the core id from global subdivision indices
## void  do_iteration();
This function is used to iterate the grid once
## void  applybcs();
This function is used to apply the boundary conditions
## void  Dirichlet();
This function is used to apply Dirichlet boundary conditions
## void  Neumann();
This function is used to apply Neumann boundary conditions
## void  Parallel();
This function is used to apply Parallel boundary conditions  
## void  create_types();
This function is used to create MPI_Datatypes for ghost cell and boundary layers of the subdomain
## void  communicate();
If the corresponding boolean in the bool communicating[4] array is true, then it starts a receive and a send for that side's ghost cell and boundary layers.
## void  set_sources();
This function is used to find the global sources that belong to the local domain and it stores them in the local_sources vector
## void  source(int  source_i, int  source_j, double  period, double  intensity, double  time);
This function is used to implement the source passed to it on the grid
