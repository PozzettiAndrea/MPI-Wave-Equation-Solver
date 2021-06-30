This is a Wave Equation solver

# Execution
**To compile, run** 

> make command
> 
To write an input file to the program, this is the correct format:

>imax jmax boundary_type nsources
>
>i j period intensity time

boundary_type is 0 for Parallel, 1 for Dirichlet and 2 for Neumann.

Note that you must provide i,j, period, intensity and time for each of the sources.

>100 100 2 2
>
>40 40 20 5 1
>
>70 70 20 5 2

In the above example, we are solving for a 100x100 grid with Neumann boundary conditions and two sources positioned at (40,40) and (70,70)

**To run it use it as follows:**

>mpirun -n <nprocessors> ./main $(cat input.txt)

**Output**

There is a flag in main called printing which I did not think useful to input manually or from input.txt. Change it if you want it to stop outputting .dat files.
.dat files are outputted to the output folder

# Post processing
To plot the .dat files in the output folder, use either
>python plot2d.py numberofimages

or

>python plot3d.py numberofimages

They will output to the images folder (make sure to create one if you don't have one!)

You have to manually input t_out

One plots in 2d and the other in 3d. An "animation.gif" file is also created in the images folder.

If you're getting strange results, remember to clear the output folder, run the solver again and plot then!

# Example plots:

![](./2DAnimation.gif)

![](./3DAnimation.gif)
