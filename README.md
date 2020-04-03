# ComplexAssembler

Python module to create macrocomplexes from chain pair structures via superimposition.

# Getting started

## Prerequisites

Before getting started, the user has to know that MODELLER has to be installed in their computer. It can be downloaded [here](https://salilab.org/modeller/download_installation.html).

## Installation

Clone the repository to your local machine and install the package.

```
git clone https://github.com/anasf97/SBI-Project

cd SBI-Project

python3 setup.py install

```
## Tutorial

After installing it, running ComplexAssembler is quite easy, you just have to execute the following command, followed by the arguments that are explained below.

```
assemble
```

### Command-line options

- **-i**: input directory with the pdb files. If it is not given, it takes the current directory as default.
- **-o**: output directory with the results. By default it creates a directory called "models".
- **-v**: show the log progression in the terminal.  
- **-e**: calculate the z-score and create a file for the energy profile of the model.
- **-opt**: optimize the macrocomplex.
- **-st**: to pass a dictionary with the stoichiometry of the chains.
- **-c**: if same chain have different names in different files, checks which chains are the same and names them the same

## Examples

In the folder examples, you will find some directories with chain pair PDB files so that you can give it a try. There you can also find  comparison of the models obtained with ComplexAssembler and the real ones.
