# Compounds Chromatograms - lcmsPlot Compound Discoverer Scripting Node

Renders one chromatogram plot per compound from a Compound Discoverer (CD) result
and adds a clickable **Plot** column to the Compounds table. Powered by
[lcmsPlot](https://github.com/computational-metabolomics/lcmsPlot).

Files in this folder (locatable via
`system.file("compound-discoverer-nodes/CompoundsChromatograms", package = "lcmsPlot")`):

- `CompoundsChromatograms.R` - the node script.
- `node.json` - the CD node definition.
- `README.md` - this file.

## 1. Install R and the required packages

Everything is installed **in the R that Compound Discoverer will run** (the same R
that `node.json`'s `ExecutablePath` points to).

1. Install **R ≥ 4.4.0** from <https://cloud.r-project.org/bin/windows/base/> (lcmsPlot requires R ≥ 4.4.0).

2. Install `lcmsPlot` and its dependencies. Start R and run:

   ```r
   install.packages("BiocManager")
   BiocManager::install("lcmsPlot")
   ```

3. Enable reading Thermo `.raw` files (the node reads the original `.raw` files
   directly):

   ```r
   BiocManager::install("rawrr")
   rawrr::installRawFileReaderDLLs()   # downloads Thermo RawFileReader (accept the licence)
   rawrr::installRawrrExe()
   ```

4. *(Development only - skip for normal use.)* To load an `lcmsPlot` **source tree**
   instead of the installed package, also install `pkgload`:

   ```r
   install.packages("pkgload")
   ```

   Then set the **lcmsPlot Source Directory** node parameter to the source path.

   **That parameter replaces the `lcmsPlot` package itself, not its
   dependencies.** Everything in the package's `Imports` still has to be
   installed in this R, so run step 2's `BiocManager::install("lcmsPlot")` at
   least once even when you intend to load from source. In particular the script
   parses `node_args.json` with **`jsonlite`** *before* it can read any
   parameter, so without `jsonlite` the node cannot see the source directory you
   set - it falls back to every default and then reports that `lcmsPlot` is not
   installed. If you hit that, check the `parameters read:` line in the log.

## 2. Register the node in Compound Discoverer

1. Copy the whole `CompoundsChromatograms` folder into Compound Discoverer's scripts
   folder (administrator rights may be required):

   ```
   C:\Program Files\Thermo\Compound Discoverer 3.5\Tools\Scripts\
   ```

   (use `Compound Discoverer 3.3` / `3.4` for those versions).

2. Open `node.json` in a text editor and check the two paths in
   `ScriptProcessorArguments`:
   - `ExecutablePath` - full path to your `Rscript.exe` (the **R ≥ 4.4.0** install),
     e.g. `C:\Program Files\R\R-4.4.1\bin\Rscript.exe`.
   - `ExecutableCommandLineArguments` - the path to `CompoundsChromatograms.R` (already
     set to the `Tools\Scripts\CompoundsChromatograms` location above); keep
     `%NODEARGS%`. That is the only placeholder Compound Discoverer substitutes -
     it becomes the path to `node_args.json`, from which the script reads the node
     parameters. Anything written after it is passed to the script verbatim.

3. In Compound Discoverer, open **Help -> License Manager** and click **Scan for Missing Features**. Confirm the dialog that new features will be available after a restart.

4. **Restart** Compound Discoverer. The node appears under **Workflow Nodes → Scripting nodes** as *Compounds Chromatograms*.

## 3. Use the node

*Compounds Chromatograms* is a **post-processing** node; drag it into the
Post-Processing Node area (it is not connected to other nodes). Its **Requested Tables and Columns** are already set in `node.json`:

```
Compounds; Compounds per File; Input Files
```

The script derives the raw file paths from the exported Input Files **File Name**
column, so the `.raw` files must be readable at their original locations.

Click the node to set its parameters (below).

## Parameters

Parameters are shown in the node's settings; CD delivers the values to the script in
`node_args.json`'s `NodeParameters`. Each has a default, so the node also runs
unconfigured.

| Parameter | Default | Purpose |
|---|---|---|
| Compounds Query | `compound_rank <= 20` | which compounds to plot (a filter over the compound table; `compound_rank` 1 = most abundant by area). Empty or `NULL` plots all. See *Checked compounds* below. |
| Output Directory | next to the `.cdResult` | directory for the PNGs |
| Plot Width | `8` | plot width (inches) |
| Plot Height | `5` | plot height (inches) |
| Plot DPI | `150` | plot resolution |
| Plot Path Mode | `absolute` | cell value form: `absolute` \| `filename` \| `relative` |
| Plot Column Name | `Plot` | new Compounds column name |
| Plot Column Position After | `Name` | place the new column after this one |
| Plot Cell Renderer GUID | `EB29D794-4F2E-4785-8B80-A24D8C0FB3E4` | CD filename cell renderer |
| lcmsPlot Source Directory | *(empty → installed package)* | set only to load an lcmsPlot source tree (development) |

### Checked compounds

The node follows Compound Discoverer's usual checked-compounds convention: **if
any compound is checked in the Compounds table, only the checked ones are
plotted; if none is checked, all of them are** (subject to `Compounds Query`).
Checking compounds in Compound Discoverer is therefore enough to pick exactly
what gets plotted, and the `Compounds Query` is ignored while any box is ticked.

To combine the two, reference `checked` in the query yourself - a query that
mentions it is used exactly as written, and is never overridden:

```
checked & compound_rank <= 5
```

The check state comes from the `Checked` column of the exported Compounds table.
Result files where no compound has ever been checked carry no such column; the
node then just reports that in its log and plots per the query. Referencing
`checked` in the query in that case is an error rather than an empty result.

## Output

- One PNG per compound in `<result>_lcmsPlot_plots/`, named `<Compounds ID>_<label>.png`.
- A `Plot` column added to the Compounds table linking each compound to its PNG.
- A diagnostic log named `lcmsPlot-cd-node.log` next to the script. If that folder is
  not writable (e.g. under `Program Files`), the log falls back to the temp directory;
  set the `LCMSPLOT_LOG` environment variable to choose a specific path.

## Standalone run

```
Rscript CompoundsChromatograms.R <node_args.json>
```

Configuration comes from the `NodeParameters` object in `node_args.json` (defaults apply when it is absent).
