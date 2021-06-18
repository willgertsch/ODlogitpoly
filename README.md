# ODlogitpoly

## Installation
1. Install PlatEMO: https://github.com/BIMK/PlatEMO
2. Add the PlatEMO directory to MATLAB's path
3. Clone this repository or download the code to your machine.
4. Add the directory to MATLAB's path.

## Usage
The best way to use this software is to open the project.m file in MATLAB and run it by clicking the "Run" button in MATLAB's user interface. After the genetic algorithm has run, a plot of the sensitivity function will be shown along with numerical values for the design. There are a number of options that can be specified by the user in the first part of the project.m file.

**Design options:**

* beta: a vector that contains nominal parameter values. If beta is length 2, then a linear predictor is used. If beta is length 3, a quadratic predictor is used.

* lower, upper: Numbers that set the lower and upper bounds for the design interval. May be Inf or -Inf.

* numpts: An integer number that sets the number of design points. You may have to run the code several times with different values to determine what the optimal number of designs points is.

**Algorithm options:**

* method: Controls the algorithm to use. Refer to the PlatEMO manual for a list of possible single-objective algorithms to use. 

* iterations: A number that controls the number of generations for the genetic algorithms. May not be accurate for other algorithms. See the maxFE variable instead.

* swarm: A number that controls the size of the swarm.

* proC: Probability of crossover. Controls how frequently solution vectors exchange traits.

* disC: Distribution index of crossover. Decreasing disC will mean more child solutions far from parents. Increasing will mean more child solutions close to parents.

* proM: Probability of mutation. Controls the frequency of random mutation of child solutions.

* disM: Distribution index of mutation. Decreasing disM will result in larger mutations from the original child solution. 

**Work in progress options:**

* opt: May add support for A, E, and G optimal designs.

The rest of the code below is functional and should not be changed unless you want to risk breaking things.
