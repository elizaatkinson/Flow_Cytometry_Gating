# Flow Cytometry Gating Scripts

These Python workflows process raw `.fcs` files from two flow-cytometry platforms and calculate fluorescence-positive population percentages after sequential gating. The repository includes instrument-specific scripts for:

- **Guava flow cytometer, gating version 1**
- **Guava flow cytometer, gating version 2**
- **Attune flow cytometer**

All scripts:

- read all `.fcs` files in the working directory
- load event-level channel data using `flowio`
- apply sequential gates for **cells** and **singlets** using ellipses in log-transformed scatter channels
- apply fluorescence thresholds to define positive events
- calculate percentages of positive events relative to the singlet population
- export a summary workbook named `gating_results.xlsx`
- save QC plots showing gates and fluorescence histograms

## Script overview

### 1. Guava script: gating version 1

This script is designed for Guava `.fcs` files using the following channel names:

- `FSC-HLin`
- `SSC-HLin`
- `FSC-ALin`
- `GRN-B-HLin`
- `RED-G-HLin`

It:

- removes events with non-positive values in required channels before log transformation
- optionally removes non-finite values after log transformation
- gates putative cells using an ellipse in $$\log_{10}(\mathrm{FSC\mbox{-}HLin})$$ vs $$\log_{10}(\mathrm{SSC\mbox{-}HLin})$$
- gates singlets using an ellipse in $$\log_{10}(\mathrm{FSC\mbox{-}ALin})$$ vs $$\log_{10}(\mathrm{FSC\mbox{-}HLin})$$
- defines GFP-positive and RFP-positive events using fluorescence thresholds on log-transformed Guava channels
- records event counts and percentages, including skipped files and failure reasons
- generates **flow-cytometry-style filled density contour plots** for the scatter-gating QC figures, with scatter fallback for sparse datasets

#### Gate settings

**Cell gate ellipse**

- center: `(3.6, 3.2)`
- width: `1.8`
- height: `1.2`
- angle: `50`

**Singlet gate ellipse**

- center: `(3.55, 3.5)`
- width: `1.7`
- height: `0.1`
- angle: `45`

**Fluorescence thresholds**

- GRN positive if $$\log_{10}(\mathrm{GRN\mbox{-}B\mbox{-}HLin}) > 2.4$$
- RED positive if $$\log_{10}(\mathrm{RED\mbox{-}G\mbox{-}HLin}) > 1.0$$

Equivalent approximate linear thresholds:

- $$10^{2.4} \approx 251$$
- $$10^{1.0} = 10$$

#### Outputs

- `gating_results.xlsx`
  - sheet: `Summary`
- optional QC PNG files per `.fcs` file:
  - `*_fsc_ssc.png`
  - `*_fsca_vs_fsch.png`
  - `*_histograms.png`

This version also:

- renames exported Guava file names into shorter safe names when possible
- logs skipped samples if required channels are missing, preprocessing removes all events, or fewer than `1000` cells remain after the cell gate
- uses high-resolution saved gating figures

---

### 2. Guava script: gating version 2

This script uses the **same Guava channels and overall workflow** as version 1, but with **different cell-gate parameters and a stricter red threshold**.

Required channels are again:

- `FSC-HLin`
- `SSC-HLin`
- `FSC-ALin`
- `GRN-B-HLin`
- `RED-G-HLin`

#### Gate settings

**Cell gate ellipse**

- center: `(3.7, 3.3)`
- width: `1.9`
- height: `1.5`
- angle: `50`

**Singlet gate ellipse**

- center: `(3.55, 3.5)`
- width: `1.7`
- height: `0.1`
- angle: `45`

**Fluorescence thresholds**

- GRN positive if $$\log_{10}(\mathrm{GRN\mbox{-}B\mbox{-}HLin}) > 2.4$$
- RED positive if $$\log_{10}(\mathrm{RED\mbox{-}G\mbox{-}HLin}) > 1.6$$

Equivalent approximate linear thresholds:

- $$10^{2.4} \approx 251$$
- $$10^{1.6} \approx 40$$

#### Key differences from Guava version 1

- larger cell gate ellipse
- red threshold increased from $$1.0$$ to $$1.6$$ on the log scale
- QC plotting may be disabled by default depending on the script setting
- if contour plotting has been applied to this version as well, scatter-gating plots are expected to be **filled density contour plots** rather than dot plots

#### Outputs

- `gating_results.xlsx`
  - sheet: `Summary`
- optional QC PNG files per `.fcs` file when plotting is enabled:
  - `*_fsc_ssc.png`
  - `*_fsca_vs_fsch.png`
  - `*_histograms.png`

Like version 1, this script records counts for loaded events, dropped events, gated cells, singlets, GRN-positive events, RED-positive events, and processing status.

---

### 3. Attune script

This script is for `.fcs` files exported from an Attune cytometer and uses a different set of scatter and fluorescence channel names.

Required channels include:

- `FSC-H`
- `SSC-H`
- `FSC-A`
- `BL1-H`
- `YL2-H`

It:

- log-transforms the Attune scatter and fluorescence channels
- gates cells using an ellipse in $$\log_{10}(\mathrm{FSC\mbox{-}H})$$ vs $$\log_{10}(\mathrm{SSC\mbox{-}H})$$
- gates singlets using an ellipse in $$\log_{10}(\mathrm{FSC\mbox{-}A})$$ vs $$\log_{10}(\mathrm{FSC\mbox{-}H})$$
- defines fluorescence-positive events by thresholding `BL1-H` and `YL2-H`
- saves plots for every processed file
- exports percentage-positive results to Excel
- can use **filled density contour plots** for scatter-gating figures, depending on the current script version

#### Gate settings

**Cell gate ellipse**

- center: `(4.7, 4.7)`
- width: `1.8`
- height: `0.6`
- angle: `40`

**Singlet gate ellipse**

- center: `(4.5, 4.5)`
- width: `1.5`
- height: `0.1`
- angle: `42`

**Fluorescence thresholds**

- BL1 positive if $$\log_{10}(\mathrm{BL1\mbox{-}H}) > 4.0$$
- YL2 positive if $$\log_{10}(\mathrm{YL2\mbox{-}H}) > 4.0$$

Equivalent linear thresholds:

- $$10^{4.0} = 10000$$
- $$10^{4.0} = 10000$$

#### Outputs

- `gating_results.xlsx`
- PNG files for each `.fcs` file:
  - `*_fsc_ssc.png`
  - `*_fsca_vs_fsch.png`
  - `*_histograms.png`

Compared with the Guava scripts, this Attune workflow exports a simpler summary table containing:

- `File Name`
- `Percentage BL1-H`
- `Percentage YL2-H`

## Summary of instrument and gate differences

| Script | Instrument | Scatter channels | Fluorescence channels | Cell gate | Singlet gate | Positive thresholds |
|---|---|---|---|---|---|---|
| Guava v1 | Guava | `FSC-HLin`, `SSC-HLin`, `FSC-ALin` | `GRN-B-HLin`, `RED-G-HLin` | ellipse centered at `(3.6, 3.2)` | ellipse centered at `(3.55, 3.5)` | GRN $$> 2.4$$, RED $$> 1.0$$ |
| Guava v2 | Guava | `FSC-HLin`, `SSC-HLin`, `FSC-ALin` | `GRN-B-HLin`, `RED-G-HLin` | ellipse centered at `(3.7, 3.3)` | ellipse centered at `(3.55, 3.5)` | GRN $$> 2.4$$, RED $$> 1.6$$ |
| Attune | Attune | `FSC-H`, `SSC-H`, `FSC-A` | `BL1-H`, `YL2-H` | ellipse centered at `(4.7, 4.7)` | ellipse centered at `(4.5, 4.5)` | BL1 $$> 4.0$$, YL2 $$> 4.0$$ |

## Input requirements

### Guava scripts

Place the Guava `.fcs` files in the same directory as the script.

Expected Guava channels:

- `FSC-HLin`
- `SSC-HLin`
- `FSC-ALin`
- `GRN-B-HLin`
- `RED-G-HLin`

### Attune script

Place the Attune `.fcs` files in the same directory as the script.

Expected Attune channels:

- `FSC-H`
- `SSC-H`
- `FSC-A`
- `BL1-H`
- `YL2-H`

## How the calculations work

For all scripts, the percentage-positive calculation is:

$$
\text{Percentage Positive} = \frac{\text{Number of gated positive events}}{\text{Number of singlet events}} \times 100
$$

The gating order is:

1. load all events from each `.fcs` file
2. remove invalid values where applicable
3. log-transform relevant scatter and fluorescence channels
4. apply a **cell gate** in scatter space
5. apply a **singlet gate** within the cell-gated population
6. threshold fluorescence channels to define positive events
7. calculate percentages relative to singlets
8. export summary results and QC figures

## QC plot style

Where updated versions are being used, the scatter-gating QC plots are rendered as **filled density contour plots**, which is consistent with common flow-cytometry presentation standards requesting **contour plots with outliers or pseudocolor plots**. Sparse datasets may fall back to light scatter plotting when contour estimation is not appropriate.

## Installation

Install the required Python packages before running the scripts:

```bash
pip install pandas numpy matplotlib flowio openpyxl
```

## Usage

Run the required script from a folder containing the relevant `.fcs` files:

```bash
python your_script_name.py
```

Each script writes its results into the current working directory.

## Typical outputs

Depending on the script, outputs may include:

- `gating_results.xlsx`
- `*_fsc_ssc.png`
- `*_fsca_vs_fsch.png`
- `*_histograms.png`

## Important notes and caveats

- The **two Guava scripts are not interchangeable**, because they use different cell-gate settings and different red fluorescence thresholds.
- The **Attune script is instrument-specific** and expects Attune channel names, not Guava channel names.
- Guava versions 1 and 2 explicitly remove non-positive values before log transformation; some Attune versions may need extra preprocessing if zeros or negative values are present.
- The Guava scripts skip files with too few gated cells after the cell gate if fewer than `1000` events remain.
- The Guava scripts produce a more detailed processing summary than the Attune script.
- Running multiple scripts in the same folder will overwrite `gating_results.xlsx` unless you rename outputs or separate workflows by directory.
- QC plots are useful for confirming that ellipse positions and histogram thresholds are still appropriate for each experiment.
- If a contour plot appears white with only a colorbar, ensure the current script version masks zero-density regions and links the colorbar to the filled contour object.

## Recommended workflow

1. Keep **Guava v1**, **Guava v2**, and **Attune** analyses in separate folders.
2. Place the correct `.fcs` files for each instrument in the matching folder.
3. Run the relevant script.
4. Inspect saved gate and histogram plots, especially after changing experimental settings.
5. Rename or archive `gating_results.xlsx` after each run to avoid overwriting earlier results.

## Suggested output naming practice

To avoid confusion, consider renaming the exported summary files after each run, for example:

- `gating_results_guava_v1.xlsx`
- `gating_results_guava_v2.xlsx`
- `gating_results_attune.xlsx`

This makes it easier to track which gating workflow was used for each dataset.
