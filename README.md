# FinalProject

Final project for the course Numerische Simulation at the University of Stuttgart in the winter term 2021/22.

## Overview

New features compared to our previous submissions:

1. inflow and outflow boundary conditions
2. custom complex geometries

## Installation and Prerequisites

- gcc
- libvtk7.1(p)
- cmake
- paraview for visualization
- Gimp for creating custom geometries
- MATLAB to transform the image to an input csv file for our program

## Usage

### Create custom complex geometry

We created arbitrary geometries as well as more complex scenerios as 8 bit grayscale images in [Gimp](https://www.gimp.org/). The pixel size of the image will later correspond to the number of cells used in the simulation.

### Transform the image to simulation input

Our program uses a .csv file to read in the geometry and import obstacle flags. A zero in the  csv table corresponds to a fluid cell, while a one means that the cell represents an obstacle.

To convert grayscale `.png` or `.jpeg` images to `.csv` files suitable for the simulation the small MATLAB script `image2inputcsv.m` is available.
Run from the [MATLAB](https://www.mathworks.de/products/matlab/index.html) command line `image2inputcsv(pathToImage, nameOfCSV)`. The algorithm will transform pixels with grayvalue < 256/2 to a fluid cell and the dark cells with a grayvalue > 256/2 to an obstacle cell.

The MATLAB script will also enforces the two-cells-criterion in the inner domain in an iterative fashion by deleting obstacle cells that do not fulfill it.

### Specify scenario

Simulation parameters for a scenario can be specified in the `<your_szenario.txt` document in `/ini`.
A path to the complex geometry's .csv file can be specified there.
Be careful that the dimensions of the image file and the dimensions specified in the scenario match.

Furthermore, settings for boundary conditions (SLIP, INFLOW, OUTFLOW) can be specified here for each domain side.

We provide you here with the following scenarious:

- padded lid driven cavity
- horizontal and vertical channels to test the INFLOW/OUTFLOW boundary conditions
- backward facing step in various resolutions
- bifurcating channel flow

### Build and run simulation

To build the simulation:

```bash
rm -rf build
mkdir build
cd build
cmake ..
make install
```

To run the simulation:

```bash
./finalproject ../ini/<your_szenario>.txt
```

### Visualization

Visualize the simulation result with [paraview](https://www.paraview.org).
Select the stack of the corresponding `.vtk` files in the `build/out` folder and click apply to visualize the simulation.

## Authors

We worked and debugged on the project 90% of the time together.

- **Author 1** - [janiswissinger](https://github.com/janiswissinger)

- **Author 2** - [kimkroener](https://github.com/kimkroener)

- **Author 3** - [magnusostertag](https://github.com/magnusostertag)
