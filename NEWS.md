# TriGer 0.99.0

* Submitted to Bioconductor.
* New `triger()` runs the full pipeline on two omics layers given as numeric
  matrices or `SummarizedExperiment` objects, matches samples by name, and
  returns feature names with each module.
* New `triger_association()` computes the signed correlation and -log10
  p-value matrices used by `greedy_3step()`.
* `greedy_peeling_XY_all()` and `greedy_3step()` take `max_blocks`, the
  maximum number of blocks returned (default 5), in place of `max_round`.
* `greedy_peeling_X_all()` gains `min_size`.
* Input matrices are checked for non-numeric, missing and infinite values.
* Tests compare the peeling functions with reference outputs shared with the
  Python package.

# TriGer 0.1.0

* Initial R implementation of the TriGer greedy peeling functions and plots.
