# Rocket_Nozzle_Design

This respository contains code to design a rocket nozzle using the method of characteristics (MoC).  I will be making a GUI for this in the near future.  For now, the user can change values in the *INPUTS* section in the main file.

## How To Run

To run the code, you'll need to download [MoC_Nozzle.m](MoC_Nozzle.m), [PM_EQUATION.m](PM_EQUATION.m), and [A_M_RELATION.m](A_M_RELATION.m) into the same directory.  Make any necessary changes in the MoC_Nozzle file, and then run it.

I'm going to be making a full derivation video in the near future, along with a detailed guide on how to run the code and what variables you can change.  CFD solutions will also be coming up shortly.

Without going into too much depth, there are a few easy variables that can be changed at the moment.  The first is the number of characteristics to use (*numChar*).  The more characteristics, the more accurate your solution will be, but it will take longer to solve and plot.  The ratio of specific heats (*g*) and specific gas constant (*R*) can also be changed.  The most important change is switching between setting a desired Mach number (*Me_Set*) and a desired area ratio (*Ae_At_Set*).  The one you don't want to use as an input needs to be set to zero!

For example, the default input is to set the area ratio to a known value of 3.  When we run the code, the derived exit Mach number from both the MoC and the quasi 1D analysis is 2.6374.  Now we can check the code by going back to the inputs, and changing the input Mach number to 2.6374, running the code, and checking the output area ratio, which does indeed equal 3, as it should.

## Solutions

Some solution times will be displayed in the command window, along with comparison to quasi-1D theory for A/A*.  The program will also plot the nozzle outline, along with contours of Mach number.  On the top half of the plotted nozzle, the characteristics mesh will also be plotted.

At the very bottom of the code, there is a section titled *OUTPUT GEO MESH FILE FOR GMSH*.  This will be used in a future YouTube video of mine, in which I will use the nozzle designed here as an input to a CFD program.  In order to ouput the mesh file (for use in the open-source GMSH), you need to change the variable *outputMesh* to 1.  It will generate a file named *Output_Mesh.geo* in the same directory that the code is in.
