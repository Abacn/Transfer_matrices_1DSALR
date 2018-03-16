# Transfer matrix method for one-dimensional short range repulsion and long range attraction particles

## Summary

- Language: Matlab (work with Matlab R2016a, R2017a)
- Publication: Hu, Yi, and Patrick Charbonneau. "Clustering and assembly dynamics of a one-dimensional microphase former." Soft Matter (2018). 
- URL: http://doi.org/10.1039/C8SM00315G


## File Description

- Pfunc_isobaric_3NN.m: Core algorithm. Transfer matrix for interaction with third nearest neighbor interactions. Isometric discretization.

- Pfunc_isobaric_3NN_sp.m: Transfer matrix implementing Simpson's rule for 
discretization.

- Pfunc_isobaric_3NN_dv.m: Transfer matrix with two-part discretization.

- findp.m: Find pressure for a given density

- findrho.m: Find density for a given pressure

- cdf.m: Cluster distribution function calculation.

- gapdf.m: Gap distribution funvction

- corlen.m: Correlation length calculation, also two-part discretization implemented.

- corlen_iso.m: Correlation length calculation, isometric discretization 
implemented. Test results are the same as corlen.m within numerical error.

- potential.m: The square-well-linear potential function