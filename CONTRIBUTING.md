# Contributing to lcmsPlot

Thank you for your interest in improving _lcmsPlot_. This guide describes how
the project is organised and what a contribution needs before it can be merged.
lcmsPlot is a [Bioconductor](https://bioconductor.org) package, so its
conventions follow the
[Bioconductor package guidelines](https://contributions.bioconductor.org/).

## Where to ask

- **Bugs and feature requests:** open an issue on
  [GitHub](https://github.com/computational-metabolomics/lcmsPlot/issues).
  For bugs, include a minimal reproducible example, the input type (e.g.
  `MsExperiment`, mzML files, a `.cdResult` file) and the output of
  `sessionInfo()`.
- **Usage questions:** ask on the
  [Bioconductor support site](https://support.bioconductor.org/) with the tag
  `lcmsPlot`.

All contributors are expected to follow the
[Bioconductor Code of Conduct](https://bioconductor.org/about/code-of-conduct/).

## Repositories and branches

| Branch / remote | Purpose |
|---|---|
| `devel` (GitHub) | Default and integration branch. **All pull requests target `devel`.** |
| `RELEASE_X_Y` (GitHub) | Tracks a Bioconductor release (e.g. `RELEASE_3_23`). Only bug fixes are backported here, by the maintainer. |
| `git.bioconductor.org` | The Bioconductor build system. Only the maintainer pushes here, from `devel` and the release branches. |

The Bioconductor and GitHub `devel` branches are kept in sync. Do not open pull
requests against the release branches; if a fix should be backported, say so in
the pull request.

## Workflow

1. **Open an issue** (or comment on an existing one) describing the change. The
   issue number is part of the branch name.
2. **Fork** the repository (external contributors) or create a branch from an
   up-to-date `devel` (members).
3. **Create a branch** following the naming convention below.
4. Make focused commits using **Conventional Commits**.
5. Open a **pull request into `devel`** once the
   [checklist](#pull-request-checklist) passes.

### Branch naming

```
<type>/<issue-number>/<short-kebab-case-description>
```

| Type | Use for |
|---|---|
| `feature` | New functionality (new layers, data sources, options) |
| `fix` | Bug fixes |
| `perf` | Performance improvements that do not change results |
| `refactor` | Code restructuring with no user-visible change |
| `docs` | Documentation, vignettes, examples |
| `chore` | Maintenance (dependencies, build, housekeeping) |
| `ci` | GitHub Actions workflows |

Examples from this repository:

```
feature/31/enable-lipidsearch-support
fix/45/cd-node-return-columns
perf/40/batch-rawrr-xic-extraction
refactor/43/use-featuredefs-rownames
chore/27/use-tibble
```

Branches use `feature/`, while the matching commit type is `feat`.

### Commit messages

Commits follow [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/):

```
<type>(<optional scope>): <imperative summary>

<optional body: what changed and why>

<optional footer(s)>
```

- **Summary line:** imperative mood ("add", not "added"), lower case, no
  trailing period, about 72 characters at most.
- **Types:** `feat`, `fix`, `perf`, `refactor`, `docs`, `test`, `ci`, `chore`.
- **Scope:** the area of the package that is affected. Scopes already in use
  include `raw-files`, `compound-discoverer`, `lipid-search`, `features`,
  `spectra`, `chrom-peak-rects`, `constructors` and `app`. Leave the scope out
  for changes that span the whole package.
- **Body:** explain *why* the change is needed when the summary does not make
  it obvious.
- **Breaking changes:** add `!` after the type/scope and a `BREAKING CHANGE:`
  footer describing the migration.
- **Issues:** reference them in a footer, e.g. `Closes #45` or `Refs #40`.

Examples from the history:

```
feat(lipid-search): add LipidSearch source and lp_lipid_search layer
fix(compound-discoverer): return only ID and plot columns from node
perf(raw-files): batch rawrr XIC extraction per file
docs: add runnable examples to all exported man pages
ci: pin Docker image to bioconductor_docker:3.23
```

```
feat(intensity-map)!: replace density argument with geom

BREAKING CHANGE: `lp_intensity_map(density = TRUE)` is now
`lp_intensity_map(geom = "density")`.

Closes #35
```

The pull request title follows the same format, because it becomes the commit
message when the pull request is squash-merged.

## Development setup

lcmsPlot requires R >= 4.4.0 and the Bioconductor devel release that matches
the `devel` branch. Install the development dependencies from the package root:

```r
install.packages(c("BiocManager", "remotes", "devtools"))
remotes::install_deps(
    dependencies = TRUE,
    repos = BiocManager::repositories())
devtools::load_all()
```

Some features need optional packages from `Suggests`:

- **ThermoFisher `.raw` files:** `rawrr`, followed by
  `rawrr::installRawFileReaderDLLs()` and `rawrr::installRawrrExe()`.
- **Precursor purity layers:** `msPurity`.
- **Interactive app:** `shiny`, `shinytoastr` and (for tests) `shinytest2`.
- **App styling:** the Shiny app ships precompiled Tailwind CSS in
  `inst/www/lcmsPlot.css`. Only if you change the app's classes or
  `tools/tailwind/input.css`, rebuild it with `tools/tailwind/build.sh`
  (Node >= 18) and commit the regenerated stylesheet.

## Coding style

Follow the
[Bioconductor coding style](https://contributions.bioconductor.org/r-code.html)
and the patterns already used in `R/`:

- Four-space indentation, no tabs, lines of 80 characters or fewer.
- `<-` for assignment; `TRUE`/`FALSE`, never `T`/`F`.
- `snake_case` for functions and arguments. User-facing plotting layers use the
  `lp_` prefix and return a function that the `+` operator applies to an
  `lcmsPlot()` object.
- Prefer `vapply()` over `sapply()`, and `seq_len()`/`seq_along()` over `1:n`.
- Call imported packages as `pkg::fun()`. Guard packages in `Suggests` with
  `requireNamespace()` and give an informative error if they are missing.
- **New files in `R/`** must be added to the `Collate:` field in `DESCRIPTION`,
  after the files they depend on.
- **New input types** should plug into the data adapters in
  `R/data-adapters.R` (`get_metadata()`, `get_detected_peaks()`,
  `get_feature_definitions()`, ...). Results from external software follow the
  `ExternalDataSource` pattern in `R/data-source*.R`.

## Documentation

- Document all functions with [roxygen2](https://roxygen2.r-lib.org/) (markdown
  enabled) and regenerate `man/` and `NAMESPACE` with `devtools::document()`.
  Never edit those files by hand.
- Every exported object must have a **runnable** `@examples` section; coverage
  is currently 100%. Use `@examplesIf requireNamespace("pkg", quietly = TRUE)`
  for code that needs a `Suggests` package. For third-party result files, write
  a minimal export to `tempdir()` as the existing `MZmineFeatureListsSource()`
  and `LipidSearchSource()` examples do, rather than using `\dontrun{}`.
- User-facing features should also be covered in the relevant vignette section
  (`vignettes/lcms_data_plotting.Rmd` or
  `vignettes/creating_publication_ready_plots.Rmd`).

## Tests

- Tests use [testthat](https://testthat.r-lib.org/) edition 3 and live in
  `tests/testthat/test-<source-file>.R`, mirroring the file in `R/` they cover
  (e.g. `R/data-source-mzmine.R` → `test-data-source-mzmine.R`). Shared
  fixtures belong in `helper-*.R`.
- Use the `faahKO` data or the files in `inst/extdata/standards-mzml.zip`, and
  `skip_if_not_installed()` for optional dependencies.
- Bug fixes should come with a test that fails without the fix.
- For `perf` changes, show that the extracted data is unchanged with the
  benchmark script:

  ```sh
  # on devel: record the reference outputs
  Rscript tools/benchmarks/bench-plotting.R --save /tmp/lcmsPlot-golden
  # on your branch: time the workloads and compare against the reference
  Rscript tools/benchmarks/bench-plotting.R --compare /tmp/lcmsPlot-golden
  ```

## Versioning and NEWS

lcmsPlot uses the
[Bioconductor version numbering scheme](https://contributions.bioconductor.org/versionnum.html)
`x.y.z`: `y` is odd on `devel` and even in a release, and the maintainer bumps
it when Bioconductor creates a release branch.

Every pull request that changes the package (code, documentation or examples)
must:

1. Bump the patch version `z` in `DESCRIPTION` (e.g. `1.1.16` → `1.1.17`).
2. Add a matching `# lcmsPlot x.y.z` entry at the top of `NEWS.md` that
   describes the change **for users**: what changed, how to use it, and a
   reference to the issue, e.g. `(#45)`. Put breaking changes under a
   `## Breaking changes` subheading with migration instructions.
3. Regenerate `codemeta.json` with `codemetar::write_codemeta()`, so the
   machine-readable metadata matches `DESCRIPTION` (version, authors,
   dependencies).

If several pull requests are open at once, rebase onto `devel` and take the
next free version number before merging.

## Pull request checklist

Before requesting a review, make sure that:

- [ ] The branch name and pull request title follow the conventions above, and
      the pull request links its issue.
- [ ] `devtools::document()` has been run and `man/`/`NAMESPACE` are committed.
- [ ] `devtools::test()` passes.
- [ ] `rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "warning")` passes
      with no new NOTEs.
- [ ] `BiocCheck::BiocCheck(".")` reports no new ERRORs or WARNINGs.
- [ ] `DESCRIPTION` version is bumped, `NEWS.md` has an entry, and
      `codemeta.json` is regenerated.
- [ ] New exported functions have runnable examples and tests.

Continuous integration (`.github/workflows/R-CMD-check.yaml`) runs the tests
and `R CMD check` in the `bioconductor/bioconductor_docker` image on every pull
request into `devel`, and must pass before merging. Vignettes are built only in
the weekly scheduled run, so build them locally
(`devtools::build_vignettes()`) if you change one.

Pull requests are reviewed by the maintainers and squash- or rebase-merged into
`devel`.
