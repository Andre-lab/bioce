## Method overview

A method based on Bayesian statistics that infers conformational ensembles from a structural library generated by all-atom Monte
 Carlo simulations. The first stage of the method involves a fast model selection approach based on variational Bayesian inference
that maximizes the model evidence of the selected ensemble. This is followed by a complete Bayesian inference of population weights in the selected ensemble.

![Alt text](images/method_pipeline.png)

## Installation
1. Download the latest 64-Bit Python 3.6 Installer from [Anaconda](http://continuum.io/downloads) and run it.
2. If you don't have git installed, follow instructions at [Git page](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git),
otherwise simply run in terminal:
```
git clone https://github.com/Andre-lab/bioce.git
```
3.	Install dependencies using yml file
```
cd bioce
conda env create -f bioce.yml
```
4. Activate conda enviroment
```
source activate bioce
```
5.	Build and install software (use –user flag if you want to install it for just single user)
```
python setup.py install
```
6. Check if scripts start up
```
python variationalBayesian.py –help
python fullBayesian.py --help
```
If you see no errors but options menu pops up, you are good to go.

## Running examples

### Simple run example
The simple example for running variational and full Bayesian approach with all necessary input files can be found in data folder
```
cd data
python ../fullBayesian.py -p flat_weights5models.txt -s TrmSimulatedIntensities5models.dat -e synthetic_60p.dat -f names5models.txt
```
After running both variational and complete Bayesian inference you should get simillar results.

### Generating input data

The above example assumed that all input data are in the right format. This may however not be the case when you start from a set of PDB models
We assume here that you have generated set of structural models and you have experimental scattering on NMR chemical shifts data available.\
You can refer to The typical workflow may look as follows.

1. Run script to generate scattering curves
This requires foxs to be installed and avaialable through PATH
```
prepareBayesian.py -s strcuture_lib_dir -e experimental_data
```
Apart from SimulatedIntensities.txt, which contains tabularized intensities for each model,
a file with starting weights and a list of files is generated. These files are needed
to run Bayesian inference in the next step.

2. Run variational Bayesian inference
Once you prepare input files as described in step 1, you can run model selection using variational Bayesian inference:
```
python ../variationalBayesian.py -p flat_weights5models.txt -s TrmSimulatedIntensities5models.dat -e synthetic_60p.dat -f names5models.txt
```


