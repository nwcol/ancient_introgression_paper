This directory contains scripts, data, and analysis for a single realization of
genome-scale coalescent simulation with a realistic demographic model and a 
heterogeneous (empirically sourced) mutation map (windowed at the 1kb scale). 

The goal is to showcase the mutation-weighting strategy: how closely can we get
data generated under a heterogenous map to conform to the pure one-locus 
expectation?

We use empirical recombination maps to achieve the same distributions of locus 
pairs across bins that we observe in real data.

the "reparse" scripts in scripts/ were written after a revision of the 
overhang parsing parameter in dpluspy and can in future reproductions of this 
analysis be ignored.

The term "weighting scheme 1" here refers to a replacement of the naive
denominator (which is the tally of two-locus pairs among considered loci) 
with a quantity which is the sum over avg(ul * ur) tabulated separately in each
interval. "Weighting scheme 2" pools ul * ur across the genome and computes 
this average once (the word used in the code for this scheme is "aggregate").
