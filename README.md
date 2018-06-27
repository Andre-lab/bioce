# BIoCE
## Bayesian inference of conformational ensembles

##Method overview
![Alt text](images/bioce_overview.png)
A method based on Bayesian statistics that infers conformational ensembles from a structural library generated by all-atom Monte
 Carlo simulations. The first stage of the method involves a fast model selection approach based on variational Bayesian inference
that maximizes the model evidence of the selected ensemble. This is followed by a complete Bayesian inference of population weights in the selected ensemble.

## Installation
1. Download the latest 64-Bit Python 3.6 Installer from [http://continuum.io/downloads] and run it.
2. If you don't have git installed, follow instructions at [https://git-scm.com/book/en/v2/Getting-Started-Installing-Git], otherwise simply run :
```
git clone https://github.com/Andre-lab/bioce.git
```
3.	Install dependencies using yml file
```
cd bioce
conda env create -f bioce.yml
```
4.	Build and install software (use –user flag if you want to install it for just single user)
```
python setup.py install
```
5. Check if scripts start up
```
python variationalBayesian.py –help
python fullBayesian.py --help
```
If you see no errors but options menu pops up, you are good to go.

## Running examples
```
cd data
python ../fullBayesian.py -p flat_weights5models.txt -s TrmSimulatedIntensities5models.dat -e synthetic_60p.dat -f names5models.txt
```
