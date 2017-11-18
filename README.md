# READ THIS!

NOTE: THIS CODE DOES NOT WORK FOR SPECIFIC HEAT RATIOS OTHER THAN 1.4 AT THE MOMENT.  The problem is in the Mach number solution in the PM_EQUATION function.  I have hard-coded g = 1.4 by accident.  This can be easily changed.  I'm currently updating the file to fix that mistake and increase flexibility for the user inputs.  Stay tuned.

# Rocket_Nozzle_Design

This respository contains code to design a rocket nozzle using the method of characteristics.  I will be making a GUI for this in the near future.  For now, the user can change values in the *INPUTS* section in the main file.

## How To Run

To run the code, you'll need to download both [MoC_Nozzle.m](MoC_Nozzle.m) and [PM_EQUATION.m](PM_EQUATION.m) into the same directory.  Make any necessary changes in the MoC_Nozzle file, and then run it.  

## Solutions

Some solution times will be displayed in the command window, along with comparison to quasi-1D theory for A/A*.  The program will also plot the nozzle outline, along with contours of Mach number.  On the top half of the plotted nozzle, the characteristics mesh will also be plotted.
