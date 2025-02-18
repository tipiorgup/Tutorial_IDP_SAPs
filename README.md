# Tutorial IDPs SAPs

## Description
This tutorial uses an adaptation of the CALVADOS2 model by 

In the folder WaNos there are several different WaNos: **DFT-Turbomole**, **DFTBplus**, **Mult-It**, **NN-Delta-ML**, **ORCA**, **Super-XYZ**, **Table-Generator** and **UnpackMol**, used to build the workflow. Below we describe each one and the main parameter exposed.

## In this tutrial, we will be able to:

1. Load a set of molecular trial structures in a `.tar` file.


## Workflow 

## 1. Installation and dependencies
Install..
### 1.1 Conda Environment and Python dependencies
Install conda environment, the following packages would be needed:

```
conda create --name environment_name python=3.6 --file environment.yml
```
Install via pip ordered-enum

```
pip3 install ordered_enum
```
## Install CALVADOS
> :warning: **If you are installing DFTB+ on a different machine**: Be very careful and make the necessary adjustments.

1. git clone -b machine-learning https://github.com/tomaskubar/dftbplus.git 
2. cd dftbplus
3. module load gnu8/8.3.0
4. module load openblas/0.3.7
5. module load cmake
6. mkdir _build 
7. FC=gfortran CC=gcc cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb+ -B _build .
8. cmake --build _build -- -j 
9. cmake --install _build

```diff 
+ Check if the `dftb+` executable exist in the dftbplus/_build/prog/dftb+/ folder. If so, then everything is okay. 
```

```diff 
+ be cautious with conda env., it must have libgfortran5 if the installation wants to be done inside the conda environment
```
## 2. Inputs and Outputs
### 2.1 Inputs
  - Molecular geometry far from the equilibrium region (as shown in **Fig 1 (a)**).
  - Set of randmized structures around the Potential energy surface (`.tar` file as shown in **Fig 1 (b)**).
### 2.2 Outputs
  - Machine learning report 
  - Pickle files of the machine learning model.

In the video **$\Delta$-learning workflow**, we teach how to set up the workflow to run it in SimStack. 

## References

1. J. Behler and M. Parrinello, Physical Review Letters 98, 146401 (2007).
2. J. Zhu, V. Q. Vuong, B. G. Sumpter, and S. Irle, MRS Communications 9, 867 (2019).
3. R. Ramakrishnan, P. O. Dral, M. Rupp, and O. A. von Lilienfeld, Journal of Chemical Theory and Computation 11, 2087 (2015).
4. C. R. Rêgo, J. Schaarschmidt, T. Schlöder, M. Penaloza-Amion, S. Bag, T. Neumann, T. Strunk, and W. Wenzel, Frontiers in Materials , 283.
5. C. L. Gómez-Flores, D. Maag, M. Kansari, V.Q. Vuong, S. Irle, F. Gräter, T. Kubař, and M. Elstner, Journal of Chemical Theory and Computation 18, 1213 (2022).
