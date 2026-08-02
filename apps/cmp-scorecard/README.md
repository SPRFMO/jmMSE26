# JM MSE candidate scorecard explorer

This Shiny application reads the current reference-OM scorecard input from
`doc/data/candidates/candidate_scorecard_input_reference.csv` and provides:

- selectable CMPs and performance metrics;
- equal, balanced, and square-root-CV dispersion weighting;
- editable custom metric weights;
- an ordered CMP score plot and results table; and
- a performance quilt recalculated for the selected CMP set.

Run from the repository root with:

```r
shiny::runApp("apps/cmp-scorecard")
```

Scores are direction-aware min--max relative preferences. They are sensitive
to the selected CMP set and are not absolute acceptability scores or agreed
management preferences.
