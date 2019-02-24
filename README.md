# GAClusterKriging

This project refers to the scripts used in the paper
"__A New Methodology for Automatic Cluster-Based Kriging using K-Nearest Neighbour and Genetic Algorithms__"

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for testing purposes.

### Prerequisites

It is advisable that the user download the RSTUDIO Desktop Version, which is FREE and OPEN-SOURCE LICENSED, available in the link:
https://www.rstudio.com/

The version 1.2.931 was used in the development of the scripts.

### Installing
The installation is simple and can be done with all default selections by pressing the Next Button.
After the installation is complete, it is necessary to install the packages included in the scripts of the project.

To install packages, select the "Tools" option in the menu and "Install Packages". 

![Alt text](/img/toolsInstallPackages.png?raw=true "Optional Title")

In the text field "Packages", type the following packages and press the button "Install". 

![Alt text](/img/toolsInstallPackages2.png?raw=true "Optional Title")

Wait for the installation to finish before installing another package.

 - gstat   
 - sp  
 - automap
 - GA
 - reshape
 - NISTunits
 - SearchTrees
 - RGeostats
 - fpc
 - outliers
 - scales

After all packages are installed, the databases must be prepared/downloaded.

1. Meuse -> Available in the "sp" package
2. Coalash -> Available in the "gstat" package
3. Broomsbarn ->  http://www.kriging.com/datasets/ (Arquivo BroomsBarn.dat) or in this repository
4. Wolfcamp -> http://www.kriging.com/datasets/ (Arquivo Wolfcamp.dat) or in this repository
5. Walkerlake -> Available in the "gstat" package

"__Broomsbarn.dat__" and "__Wolfcamp.dat__" must be included in the same path/directory of the script (GAClusterKriging.R) file

After theses steps, you may run the code!

## Running the tests

Open the script (GAClusterKriging.R) File in the RSTUDIO software.

__BEFORE RUNNING ANY TESTS, MAKE SURE TO FOLLOW THESE STEPS:__

 - In the toolbar select "Session", "Set Working Directory", "To Source File Location";
 - In the toolbar select "Session", "Clear Workspace";

To RUN the tests, first the database must be selected in the section bellow!

![Alt text](/img/code1.png?raw=true "Optional Title")

Next, type the name of the output file and the number of tests that will be applied to each cluster.

![Alt text](/img/code2.png?raw=true "Optional Title")

Next, select the size of the Genetic Algorithm Population and Number of Iterations.
Next, select the number of neighbours considered in the KNN method.
Is important to note that the number is N-1, because the closest point, is the point itself.
So if you want 3 neighbours, you must type 4.
Finally, select the number of clusters of the K-Means Method.

![Alt text](/img/code3.png?raw=true "Optional Title")

Finally, save the file and click on the "Source" Button to run the script.

![Alt text](/img/code4.png?raw=true "Optional Title")

The execution can be time demanding depending on the database.

### Tests Results

The output text file will be created in the directory of the R script file. 

At the header, informations about the selected parameters are shown.

Then, iterations for each number of clusters are listed. 

![Alt text](/img/outputfile.png?raw=true "Optional Title")


## Authors

* **Carlos Yasojima**

takeshiyasojima@gmail.com
