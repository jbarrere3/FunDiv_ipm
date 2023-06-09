# FunDiv_ipm

Repository containing the code to run the model and statistical analyses of the paper "_The functional composition of tree communities interacts with climate to shape forest resilience to storm disturbances_" by Julien Barrere, Björn Reineking, Maxime Jeaunatre and Georges Kunstler, currently in preparation.

The code requires prior installation of the [matreex](https://github.com/gowachin/matreex) R package, developped by Maxime Jeaunatre (INRAE), and of the ```targets``` package. The data folder, required to run the code, can be made available upon request to julien.barrere@inrae.fr

Once the packages are installed and the data folder is placer in the main folder, just run ```targets::tar_make()``` from R and the script will download the other packages required and run the analyses. 

