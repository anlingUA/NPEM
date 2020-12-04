This function peforms nonparametric mediation analysis by capturing unusually high contributed
information scores for each mediation pathway. All feature sets must be in a matrix or data frame
format.

Currently, four algorithms are available:
1) UV  - iterative selection for continuous, discrete, or ordered mediators with significance
	 estimated at each iteration

2) UVS - a modified version of UV with significance estimated after all iteration steps

3) BV  - iterative selection for hurdle model mediators with a binary hurdle for presence-absence
	 with significance estimated at each iteration

4) BVS - a modified version of BV with significance estiamted after all iteration steps
