# Biotic Novelty Index: index calculation and simulation

This repository provides the base code and example applications for the following article: 

  **Conrad Schittko** & **Maud Bernard-Verdier**, Tina Heger, Sascha Buchholz, Ingo  Kowarik, Moritz von der Lippe, Birgit Seitz, Jasmin Joshi  & Jonathan M. Jeschke.  A multidimensional framework for measuring biotic novelty: How novel is a community? (2020) *Global Change Biology*, accepted on April 15th 2020 (first submitted December 2019)

The code makes use of R packages "vegan" (version 2.5-5) and "FD"(version 1.0-12), and dependencies.

## Content of the repository:

- *FUNCTION BNI.calc.R* : function calculating the Biotic Novelty Index (BNI)

- *FUNCTION BNI simulation functions.R* : a set of functions to simulate communities and explore the behaviour of the BNI, and two functions for plotting results 

- *Simulation testing.R* : R script to re-run the 9 different simulation scenarios from the article

- *Case study 1.Rmd* : R markdown document replicating the case study 1 on anonymized data

- **/data** : example data folder
    - *berlin.rda* : anonymized data for case study 1