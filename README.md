# Quarto Document Template

This repository is a minimal Quarto workspace with:

- `index.qmd`: starter document wired to `references.bib` and the AMA CSL, including demo paragraphs that cite `@smith2023`, `@lee2022`, `@chen2021`, and `@garcia2020`.
- `references.bib`: sample BibTeX database you can replace with your own sources.
- `AMA.csl`: American Medical Association CSL style downloaded from Zotero.
- `helper.R`: helper script that validates citations in `.qmd` files against the BibTeX file and optionally checks DOI resolvability.

## Getting Started

1. Install Quarto (https://quarto.org/docs/get-started/) and make sure `Rscript` is available on your path.
2. Clone or copy this folder to wherever you plan to write your document.
3. Edit `index.qmd` (or add more `.qmd` files) with your content and cite using `@citekey` syntax.
4. Replace `references.bib` with your bibliography. You can paste BibTeX entries directly into the file, one per reference.

## Demo Content

- The Introduction shows inline citations and grouped citations so you can see how AMA formatting looks.
- The Workflow Demo section walks through the helper script + render process.
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
