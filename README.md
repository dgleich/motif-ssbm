Numerical Analysis of Motif-weighted adjacency matrices for SSBM
================================================================

This repository reproduces some of the new experiments shown at
the Biennial Numerical Analysis Conference in my (@dgleich) 
presentation on higher-order analysis of networks based on 
motifs. 

I've added all the codes, but the key ones are as follows.

* `motif_codes.jl` basic functions to generate SSBMs and
some simple analysis
* `matrixops.jl` ways to generate standard matrices based
on an adjacency matrix (or weighted adjacency)
* `info-threshold.jl` compute the various information 
theoretic thresholds in terms of the SSBM (based on Abbe
et al., community detection in the stochastic block model)
* 

Helpful scripts to reproduce various figures.

## Illustrative introduction figures
* `motif-weighting-pictures.jl` make images of the SSBM adjacency matrices
* `gnr_figure.jl` make some nice looking pictures of graphs in the plane
* `apache_figure.jl` show a nice picture of a mesh to illustrate sci-comp view
* `animate_ssbm.jl` produce an animation of a single SSBM generation as the graph fills in

## Analysis figures
* `powermethod_steps.jl` Make the picture showing that the power
method does better with motif-weighting at recovery
* `powermethod_animation.jl` Figures of what happens as you run the
power method 
* `eigenvalue_sensitivity.jl` look for a reason why the powermethod
converges faster for motifs.
* `eigenvalue_histograms.jl` show histograms for where the spectrum lies
* `eigenvector_sharpness.jl` show that the eigenvectors for the motifs are
sharper, which illustrates our conjecture about the second reason that results
are better for the motif-weighting

## Working codes
* `acl_animation.jl` and `acl_steps.jl` shows what happens with the ACL
algorithm instead of the power method. 
* `expected_eigenvalues.jl` shows some guesses at analytic forms for 
the eigenvalue regions and gaps
* `powermethod_steps_clique4.jl` shows the results using 4-clique weighted
matrices instead of triangles. 
* `Marchenko-Pastur-law-example.jl` testing for a Marchenko-Pastur law for the motif-weighting
* `Beyond-Marchenko-Pastur-law-example.jl` trying to guess a probability law for the remainder