# FinalProject

Final project for the course Numerische Simulation at the University of Stuttgart in the winter term 2021/22.

## Content

1. complex geometries
2. boundary conditions
3. ?

## Installation and Prerequisites

- gimp for image creation
- matlab for image preprocessing
- gcc
- libvtk7.1(p)
- cmake
- paraview for visualization

## Usage

### Create custom complex geometry

To create a black-white image of your custom complex geometry e.g. in [Gimp](https://www.gimp.org/).
Notice that the image has to respect the *two-cells-criterion*: very obstacle has to be at least two cells thick in x- and in y-direction.

To convert black-white `.png`or `.jpeg` images to `.csv` files suitable for the simulation:

Run from the [matlab](https://www.mathworks.de/products/matlab/index.html) command line: `imageToInputCSV.m <your_image>`

The `matlab` script will enforce the two-cells-criterion.

### Specify scenario

The scenario can be specified in `<your_szenario.txt` document in `/ini`.
Further, the path to the complex geometry can be specified there.
Be careful that the dimensions of the image file and the dimensions specified in the scenario match.

An exemplary scenario is preconfigured in ?.

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
Select the stack of the corresponding `.vtk` files in the `build/out` folder and click apply.
Then you can play the simulation.

## Authors

We worked on the project mostly together.

* **Author 1** - [janiswissinger](https://github.com/janiswissinger)

* **Author 2** - [kimkroener](https://github.com/kimkroener)

* **Author 3** - [magnusostertag](https://github.com/magnusostertag)
