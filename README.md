# jmMSE26

Documentation and implementation guidance for candidate management procedures
(CMPs) evaluated for SPRFMO jack mackerel.

This repository is a reporting and specification layer derived primarily from
the current `jmMSE` implementation and outputs. Scientific model development
and MSE execution remain in
[`SPRFMO/jmMSE`](https://github.com/SPRFMO/jmMSE). The site reorganizes the
current candidate definitions, workflow, test evidence, and implementation
instructions. SCW17 material is retained as historical context.

## Render the site

```sh
cd doc
quarto render
```

Rendered files are written to `docs/` for GitHub Pages.

## Repository roles

- `doc/cmp/`: CMP registry, component inventory, and specification template.
- `doc/application/`: annual data-update and application procedures.
- `doc/workflow/`: relationship between assessment, MSE, and management advice.
- `doc/evidence/`: evidence supporting CMP screening and selection.
- `doc/governance/`: provenance, versioning, and change-control rules.
- `doc/data/cmp-registry.csv`: machine-readable registry sourced from `jmMSE`.
- `R/build_candidate_slick.R`: validated eight-CMP performance-to-Slick adapter.
- `output/jm_candidates.slick`: the single validated eight-CMP review file,
  combining the om11 reference with the single-stock and two-stock robustness
  sets (100 simulations, 1970–2050), with Set, OM, and Stock filters and
  one-click OM presets.

Candidate status is explicit. A test candidate is not an adopted management
procedure unless its registry status and decision record say so.

## Open the current Slick results

The validated file is
[`output/jm_candidates.slick`](output/jm_candidates.slick). It contains eight
CMP comparison cases across the reference and robustness operating models.
Northern and Southern stock components retain their source labels. For readable
plots, F/FMSY is missing where FMSY is zero and finite F/FMSY values above 4
are displayed at 4; the original `jmMSE` performance tables are unchanged.

### Use the hosted Slick app

1. [Download `jm_candidates.slick`](https://raw.githubusercontent.com/SPRFMO/jmMSE26/main/output/jm_candidates.slick).
2. Open the [Blue Matter Slick app](https://shiny.bluematterscience.com/app/slick).
3. If the service says the app has stopped, select **Restart app**.
4. On the opening page, under **Load your MSE Results**, select the downloaded
   `.slick` file.
5. Wait for the Overview page, then use the MP and OM filters to explore the
   results.

The hosted app accepts a file upload but does not currently provide a field for
loading a Slick file directly from a GitHub URL.

### Open it from R

Install `Slick` once if needed, then download, check, and open the file:

```r
install.packages("Slick")
library(Slick)

url <- paste0(
  "https://raw.githubusercontent.com/SPRFMO/jmMSE26/",
  "main/output/jm_candidates.slick"
)
file <- tempfile(fileext = ".slick")
download.file(url, file, mode = "wb")

slick <- readRDS(file)
Check(slick)
App(slick = slick)
```

### Add it to the hosted app's case-study menu

The hosted app obtains its selectable case studies from the
[`Blue-Matter/SlickLibrary`](https://github.com/Blue-Matter/SlickLibrary)
repository. To request permanent menu integration:

1. fork or clone `Blue-Matter/SlickLibrary`;
2. copy this file to
   `Slick_Objects/SPRFMO_Jack_Mackerel.slick` in that repository;
3. submit a pull request to Blue Matter; and
4. after it is merged, select **SPRFMO Jack Mackerel** from the app's case-study
   menu.

This last step requires approval from the SlickLibrary maintainers. The
[official Slick developer guide](https://slick.bluematterscience.com/articles/DevelopersGuide.html#slick-object-library)
describes the same contribution process.
