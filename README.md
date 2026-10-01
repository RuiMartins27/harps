# HARPS

HARPS (Harmonic Absorbed Radiation Plasma Solver) is a C++ module for calculating the absorbed power density profile and the electric field from a microwave source inside a plasma reactor.

The goal of the code is to couple this module with a plasma model, that being a reduced dimension or CFD model.

## Installation

The code uses the PETSc library to solve the linear complex system of equations that arrives from the fields. This package is obtained from the spack package manager. The installation guide for this package is the following.

```bash
git clone https://github.com/spack/spack.git ~/spack
source ~/spack/share/spack/setup-env.sh
spack install petsc +complex clanguage=C++
```

Besides the installation step, the enviromental variables need to be set, either by running the follwoing commands everytime the machine is restarted or adding these commands to the .bashrc file

```bash
source ~/spack/share/spack/setup-env.sh
spack load petsc
export LD_LIBRARY_PATH=$(spack location -i petsc)/lib:$LD_LIBRARY_PATH
```

## How to run

The current compilation is done simply using a Makefile and the code is ran using a single command line.

```bash
make
```

Here the input file is the example1.in which is a waveguide WR-340 with a random source which results on the propagation of the only mode that can propagate inside this waveguide (TE10)). There is also the template.in that includes all possible input options. The other included examples are antenna on a uniform plasma, wave propagating and encounting a plasma besides others.


```bash
mpirun -np 1 ./harps.exe input/example1.in
```

The analysis is currently python codes. The most important one is to graph the resulting fields


```bash
python3 analysis/3D_graph.py
```

or with the added 'y' input for 2D slices in the y direction instead of the z direction

```bash
python3 analysis/3D_graph.py y
```

## Authors

This project is being developed by the Circular Chemical Engeneering department of Maastricht University. This project is personally being developed mainly by Rui Martins (rui.piresmartins@maastrichtuniversity.nl).


## How to cite

Harps has been described in a recent publication and can be cited using Rui Martins et al 2026 Plasma Sources Sci. Technol. 35 095012


## Contributing

If you find a bug please contact the authors.

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## When coupling with PeleLMeX installation guide is as follows

Install Dependencies

```
git clone https://github.com/spack/spack.git ~/spack
source ~/spack/share/spack/setup-env.sh
spack install petsc +complex +mumps +metis +scalapack +hypre +suite-sparse +superlu-dist +fortran +shared ^mumps+parmetis
spack install openmpi
```

Clone the directories
```
mkdir power_coupling
cd power_coupling
git clone git@gitlab.com:Ruimartins27/harps.git
git clone --branch axiswirl-power --single-branch https://github.com/RuiMartins27/PeleLMeX.git
```


Clone the PeleLMeX Submodules
```
cd PeleLMeX
git submodule update --init --recursive
```


Install Chemistry Set if needed. After adding new chemistry set to PelePhysics/Mechanisms such as airThermal.
When installing airThermal it is common to get error on line 3285 of airThermal/mechanism.H where it is missing amrex::Real  before redP
```
poetry install
poetry update
cd Submodules/elePhysics/Support/ceptr/
poetry run convert -f ../../Mechanisms/airThermal/mechanism.yaml
```


Set Environmental Variables
```
source ~/spack/share/spack/setup-env.sh
spack load petsc
spack load openmpi
export LD_LIBRARY_PATH=$(spack location -i petsc)/lib:$LD_LIBRARY_PATH
```

Compile
```
make TPL
make -j4
```

The input file used by harps when called inside PeleLMeX is harps/input/2D_RZ.in
Other input files are created using harps/analysis/write_input_files.py . Namely the problem == 15 is the one used to create harps/input/2D_RZ.in (although it is already pre-created in git).  Other useful analysis scripts can be found inside analysis directory.