Simulating simple H-M community selection, selecting top 2, fixed total Newborn biomass.
* In Heritability/ , heritability of community function and determinants are calculated from C549 of the simulation.
* With Heritability/NoSpike/DeterminantDist.m , fp(0) and phiM(0) of all Newborns are normalized to compare to the normal distribution.
  * phiM(0) distribution is normal, but not fp(0) distribution. For example, parent adults(40) has 5 null mutant while the other fp values are quite uniform. Offspring newborns(2) inherits 2 null mutants, its fp(0) = 0.1393; Offspring newborns(11, 27, 56) inherits 1 null mutant, their fp(0) ~ 0.143. These 4 offspring deviates significantly away from the mean, resulting in the heavy tail in the negative region.
