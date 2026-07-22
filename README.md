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
- `R/build_candidate_slick.R`: validated MP29/MP32 performance-to-Slick adapter.

Candidate status is explicit. A test candidate is not an adopted management
procedure unless its registry status and decision record say so.
