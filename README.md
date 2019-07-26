## Occupancy modeling with RStan

This package uses RStan for occupancy modeling.

The package is currently being adapted from a [tutorial](https://github.com/rerickson-usgs/StanOccupancyModelTutorials/). 

## Project status and to do:

Currently, several Stan models exist. 
Most wrapper functions have not yet been written. 


- [   ] Finish writing R wrapper, including placing more code inside R functions and writing examples 
 - [   ] `USFWS_UMR_eDNA.R` (need to move code inside function)
 - [   ] `occupancy.stan`
 - [   ] `eDNAoccupancy.stan`
 - [   ] `occupancyBernoulli.stan`
 - [   ] `occupancy_bernoulli.stan`
 - [ x ] `lm.stan`
 - [   ] `occupancy_mu.stan`
 - [ x ] `logistic_stan.stan`
 - [   ] `occupancy_mu_binomial.stan`
 - [ x ] `logistic_target_stan.stan`
- [   ] Write test functions (while doing above)
- [   ] Create plotting functions 
- [   ] Upload example datasets 
- [   ] Move tutorial into Vignette
- [   ] Update this file 

My plan of attack:

1. Start working through tutorial files and moveing them over here.
2. Write unit testing while doing this.
3. Do this in order because this will be easiest.
4. Convert existing tutorials into vignettes for this project. 

## Contact

Richard A. Erickson (rerickson@usgs.gov)

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

This software is provided "AS IS".
