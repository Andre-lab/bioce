# bioce
bayesian inference of conformational ensembles

Place holder for bayesian inference method for problems like two domains connected by flexible linker 

Conda enviroment can be setup with:
```
conda env create -f bioce.yml
source activate bioce
```
And to run the script
```
cd data
python ../fullBayesian.py -p flat_weights5models.txt -s TrmSimulatedIntensities5models.dat -e synthetic_60p.dat -f names5models.txt
```
