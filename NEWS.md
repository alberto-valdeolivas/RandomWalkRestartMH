## RandomWalkRestartMH v1.13.2 (Release date: 2026-08-04)

* Removed the archived `dnet` dependency. Top-result network construction now
  uses an internal `igraph` implementation with equivalent behavior for the
  subset of `dNetInduce` functionality used by the package.

* Fixed inconsistent multiplex node-count field access in heterogeneous network
  construction and multiplex random-walk validation.

* Improved validation of transition and restart parameters.

* Fixed custom layer weights (`tau`, `tau1`, and `tau2`) to use the same
  sum-to-one convention as their defaults, and fixed validation of seed nodes
  from the second multiplex network.

* Updated deprecated `igraph` calls used by package code for compatibility with
  current `igraph` releases.

## RandomWalkRestartMH v1.9.2 (Release date: 2020-09-01) 

Minor update to fix some problems with the *create.multiplexHetNetwork.topResults*
function when using weights in the bipartite network. 

## RandomWalkRestartMH v1.9.1 (Release date: 2020-08-05) 

Major update with the following main points:

* Added the possiblity to perform Random Walk with restart on full multiplex
heterogeneous networks (Both networks connected by bipartite interactions
can now be multiplex).

* RandomWalkRestartMH can now handle weighted networks. 

* New package website for the sake of clarity. 

* Minor changes in the way to create multiplex objects (A list of graphs 
is now required.)

* Replaced the references to bioRxiv by the final publication in Bioinformatics.
