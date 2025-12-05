# Quarto Document Template

## About

This repository packages a minimal Quarto project plus an R-based reference validator so you can prototype manuscripts, verify citations/DOIs, and publish directly to GitHub Pages. Use it as a starting point for new reports or as a sandbox for testing Quarto workflows.

This repository is a minimal Quarto workspace with:

- `index.qmd`: starter document wired to `references.bib` and the AMA CSL, including demo paragraphs that cite `@smith2023`, `@lee2022`, `@chen2021`, and `@garcia2020`.
- `sections/`: modular `.qmd` files (Abstract, Introduction, Methods, Results, Discussion, Conclusion, Acknowledgements) that `index.qmd` imports via Quarto include shortcodes.
- `references.bib`: sample BibTeX database you can replace with your own sources.
- `AMA.csl`: American Medical Association CSL style downloaded from Zotero.
- `helper.R`: helper script that validates citations in `.qmd` files against the BibTeX file and optionally checks DOI resolvability.

## Getting Started

1. Install Quarto (https://quarto.org/docs/get-started/) and make sure `Rscript` is available on your path.
2. Clone or copy this folder to wherever you plan to write your document.
3. Edit `index.qmd` (or add more `.qmd` files) with your content and cite using `@citekey` syntax.
4. Replace `references.bib` with your bibliography. You can paste BibTeX entries directly into the file, one per reference.
5. Edit the files under `sections/` to update each manuscript section; `index.qmd` automatically aggregates them in IMRaD order.

## Demo Content

- The Introduction demonstrates a full CARS-style narrative with inline and grouped citations.
- The Methods, Results, and Discussion modules provide structured bullet-point prompts you can replace with study-specific prose.
- Each section currently contains placeholders for key data (design choices, numerical results, interpretations); swap them out with your manuscript content when ready.
- `references.bib` ships with four sample entries that correspond to the cite keys used in `index.qmd`.
- Use this setup to test the helper script before swapping in your own text and bibliography.

## Validate Citations and DOI Status

The helper script can scan your `.qmd` files and verify that every cited key exists in the BibTeX file. It can also issue HEAD requests to `https://doi.org/<doi>` to confirm DOI responses.

```bash
# Validate everything automatically (searches for *.qmd and references.bib)
Rscript helper.R

# Specify files explicitly
Rscript helper.R paper.qmd references.bib

# From an R session
source("helper.R")
check_refs(check_doi = TRUE)
```

The script prints:

- Missing citations (`@key` used in `.qmd` but not found in `.bib`).
- Unused bibliography entries.
- DOI validation summary (valid, invalid, missing).

## Render the Document

After the references pass validation, render the Quarto project (the front matter already requests both HTML and DOCX outputs):

```bash
quarto render index.qmd
```

Quarto writes `index.html` and `index.docx` into `_site/` by default. Use `quarto render index.qmd --to html --to docx` if you explicitly want to control the formats on the command line.

## Workflow Summary

1. Paste or write text into `.qmd` files.
2. Paste corresponding BibTeX entries into `references.bib`.
3. Run `Rscript helper.R` (or use `check_refs()` in R) to ensure citations/DOIs are valid.
4. Run `quarto render index.qmd` to generate HTML + DOCX outputs under `_site/`.

## GitHub Pages Deployment

- The workflow at `.github/workflows/quarto-gh-pages.yml` runs on every push to `main`. It checks out the repo, installs Quarto + R, installs `httr`/`bib2df`, renders `index.qmd` into `_site/`, and uploads that folder as the Pages artifact.
- In the GitHub repo settings, enable Pages with **Source → GitHub Actions** so that `actions/deploy-pages` can publish from the workflow.
- After the workflow succeeds, the `deploy` job publishes the generated site to the GitHub Pages environment and surfaces the live URL in the Actions log (https://htlin222.github.io/quarto-doc/).
- You can always browse the published document at https://htlin222.github.io/quarto-doc/.

## Helper Script Quick Reference

- Validate everything automatically (search for all `.qmd` files + `references.bib`):
  ```bash
  Rscript helper.R
  ```
- Validate specific files:
  ```bash
  Rscript helper.R paper.qmd references.bib
  ```
- From an interactive R session:
  ```r
  source("helper.R")
  check_refs(check_doi = TRUE)  # set FALSE to skip DOI checks
  ```
