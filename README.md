# jmMSE26

Documentation and implementation guidance for candidate management procedures
(CMPs) evaluated for SPRFMO jack mackerel.

This repository is a reporting and specification layer. Scientific model
development and MSE execution remain in
[`SPRFMO/jmMSE`](https://github.com/SPRFMO/jmMSE). The site reorganizes the
management-procedure inventory, workflow, test evidence, and implementation
instructions developed for SCW17.

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
- `doc/data/cmp-registry.csv`: machine-readable registry seeded from SCW17.

Candidate status is explicit. A test candidate is not an adopted management
procedure unless its registry status and decision record say so.

