# Finite Theory bench experiments: records and analysis scripts

Release **bench-v1.1** (1 November 2026). These are the records and scripts behind Section 5 (Laboratory experiments) of
*Finite Theory: Revisiting the Variable-Speed-of-Light Model of Gravitation, with Exploratory Bench Tests*, P. Bouchard,
Fornux Inc. Every number in Tables 3 to 6,
the fits of Section 5.4.2 and Figures 9, 10, 12 and 13 of the paper is reproduced by the commands below. Each folder holds
an `expected_output.txt` with the output of this release.

The release is the file `bench-v1.1.zip` on the master branch of https://github.com/philippeb8/finite-theory; the paper
gives its SHA-256 checksum.

## Contents

| Folder | Experiment | Records | Scripts | Paper |
|---|---|---|---|---|
| `freefall/` | 40 paired drops, NdFeB and steel spheres | `data/day7..day11/FreeFall-7pm-{Up,Down}.py` (tracked positions, 240 fps, in their `DATA` blocks) | `ff_refit.py` | Sec. 5.1, Table 3, Fig. 9 |
| `pendulum/` | Magnets and steel sphere on the same nylon twine | `tracks/IMG_8754_steel.csv`, `tracks/IMG_8755_magnets.csv` | `track_pendulum.py`, `periods.py`, `corrections.py` | Sec. 5.2, Table 4, Fig. 10 |
| `floating_magnet/` | Float in the afternoon, at night and on a plastic support; phone magnetometer | `tracks/*.csv`, `magnetometer/*.csv` | `track_night.py`, `track_waterline.py`, `fit_drift.py`, `magnetometer/holds.py` | Sec. 5.3, Table 5, Figs. 11 and 12 |
| `spinning_disc/` | Disc on a beam balance | readings in `readings.py` | `readings.py`, `fit_law.py` | Sec. 5.4, Table 6, Eq. (100), Fig. 13 |
| `figures/` | | | `make_figures.py` (writes `figures/out/`) | Figs. 9, 10, 12, 13 |
| `galaxies/` | Rotation-curve survey: 497 curves of SPARC and the Sofue atlas; the Milky Way | survey: read from github.com/philippeb8/grc at commit 9f3c865; Milky Way: `milky_way_sofue2020.csv` | `galaxy_fit.py`, `baseline_arctan.py`, `summary.py`, `milky_way.py`; fitted parameters in `galaxies.csv` and `milky_way.json` | Secs. 2.2.1 and 4.3.1, Eq. (91) |
| `archive/` | Interferometer experiment proposed earlier for the ISS, with its open steps | | `ISS_proposal_archived.tex` and `.pdf` | removed from the paper; a historical proposal, not a prediction of the present paper |

The gallium-vortex run (Sec. 5.5) is a single reading, 2.6 +- 3.7 g, and has no file.

## Run

Python 3.10 or later with numpy, scipy and matplotlib (`pip install -r requirements.txt`). From this folder:

```
python3 freefall/ff_refit.py freefall/data       # accelerations by half of the drop, Table 3; writes freefall/lower_half.json
python3 pendulum/periods.py                       # periods of both runs; writes pendulum/periods.json
python3 pendulum/corrections.py                   # Table 4 and the conditional readings; writes pendulum/corrections.json
python3 floating_magnet/fit_drift.py              # Table 5
python3 floating_magnet/magnetometer/holds.py     # field at each 10 cm hold
python3 spinning_disc/readings.py                 # Table 6 readings and calibration
python3 spinning_disc/fit_law.py                  # fits of Eq. (100); writes spinning_disc/fit_law.json
python3 figures/make_figures.py                   # Figs. 9, 10, 12, 13 (needs the json files above)
python3 galaxies/milky_way.py                     # Milky Way fit of Eq. (91); writes galaxies/milky_way.json
```

The galaxy survey needs the data of the GRC repository and takes a few minutes:

```
git clone https://github.com/philippeb8/grc && git -C grc checkout 9f3c8658fd7be9aa43964caf29e8ed39a0d4f575
cd galaxies && python3 galaxy_fit.py ../grc && python3 baseline_arctan.py ../grc && python3 summary.py
```

Re-tracking from the videos needs ffmpeg: `track_pendulum.py <video> <out.csv>`, `track_night.py <video> <name>`,
`track_waterline.py <IMG_8765 video> <IMG_8766 video>`. The videos are
available from the author on request; `SHA256SUMS_videos.txt` identifies the files the tracks were made from. A hash
identifies a file but does not provide its contents or validate the tracking: the tracks supplied here are the record,
and re-tracking from the videos is a separate, independent check.

## Notes on the records

- **Free fall.** `ff_refit.py` reads the `DATA` blocks of the ten session scripts without executing them. The per-session
  summary that those scripts print fits the median of each block with the release height fixed at zero; it is biased and
  is not used. Each drop is calibrated on its own steel sphere; the drop with a failed steel track (day 8, D) is excluded.
  `per_drop.csv` lists, for every drop, the calibration scale and curvature, the residuals of the calibration and of the
  two-half fit, the mid-height time and the fitted accelerations. The calibration fitted to the steel sphere is applied to
  the magnet 57 cm away in the image; the uncertainty of that transfer is not propagated, so the quoted errors describe
  the scatter of the drops (pooled, +-0.41 m/s^2) or of the five days (+-0.51 m/s^2), not a total uncertainty.
- **Pendulum.** Both records are used over the same swings: the 21 turning points at the near end and the 20 at the far
  end that follow the release. The near-end and far-end estimates come from the same recording and may share systematic
  errors, and with one recording per bob the reproducibility between runs is not measured. The tracker renews its
  background every second, which follows a shift of the camera of about 2 pixels between 52.6 and 54.5 s of IMG_8754.
  The rope and centre corrections are calculated from the geometry; the correction for the slower decay of the magnets'
  swing is measured from the peak-to-peak height of the bob in the image, scaled to the release amplitude; the stretch of
  the nylon twine is a sensitivity interval for moduli of 0.5 to 2 GPa. `corrections.py` gives the remaining difference
  two conditional readings: a gravitational-to-inertial mass ratio for equal effective lengths, or an effective-length
  difference for an unchanged mass ratio. Neither is a measured property of the apparatus.
- **Floating magnet.** The afternoon tracks (IMG_8733, IMG_8735) were made with an earlier version of the
  background-subtraction tracker, which was not kept; the tracks are the record. The six drifts and IMG_8765 are fitted
  by the same script, and their initial accelerations are outputs of its assumed constant-force and linear-drag model;
  IMG_8766, which slows down after its turn, is described by a coast (initial speed and time constant). IMG_8765 and
  IMG_8766 stood on a plastic box, with no steel under the container; their tracks follow an edge of the float in a band
  at the waterline (`track_waterline.py`).
- **Spinning disc.** The pointer positions were read from the videos; `readings.py` holds them with the coin calibration.
  `fit_law.py` makes the principal descriptive fit on the four readings of runs IMG_8693 to IMG_8695 between 2369 and
  2503 rpm, clear of the beam stop, and two sensitivity fits that add the 2591 rpm reading (within its uncertainty of
  the stop) and the provisional 2692 and 2751 rpm readings (transferred scale, stop not known). The readings were
  selected after inspecting the data; parameter errors are scaled by sqrt(chi2/dof), which does not correct clipping,
  an unknown zero or a shared calibration bias, and the uncertainty of the speeds is not included.
- **Galaxies.** `galaxy_fit.py` follows `src/main.cpp` of the GRC repository: lambda_g on a grid from 1e18 to 1e22 kg/m, m between
  1/15 and 15 times a photometric mass scale, omega (rad/s, the angular velocity of the comoving framework of the galaxy)
  set so the curve passes through the last point, objective = standard
  deviation of the residuals about their mean (unweighted). `baseline_arctan.py` fits an empirical arctangent curve with the
  same linear term, objective and parameter count. Duplicate objects are kept, as in the original survey, and listed.
  `galaxy_fit.py` uses G = 6.67408e-11, the value of the paper; `src/main.cpp` and bench-v1.0 used 6.671e-11.
  `milky_way.py` fits the Milky Way, which is not one of the 497 curves, by the same procedure: the curve of Sofue
  (2020, Galaxies 8, 37, Appendix A, Tables 5 and 6) in `milky_way_sofue2020.csv`, with the stellar mass of the Galaxy,
  5e10 Msun (Bland-Hawthorn & Gerhard 2016), as the mass scale. The principal fit starts at 2.5 kpc, outside the bulge;
  this inner limit was chosen after inspecting the fits, and the fit to the whole curve and the values for inner
  limits of 2 to 5 kpc are printed with it.

## Changes from bench-v1.0

Every file that differs from bench-v1.0 is listed here.

- **Link.** The release is the file `bench-v1.1.zip` on the master branch; bench-v1.0 pointed to a `main/bench` path.
  This archive supersedes bench-v1.0 and its unpacked `bench/` folder.
- **Pendulum** (`track_pendulum.py`, `tracks/IMG_8754_steel.csv`, `tracks/IMG_8755_magnets.csv`, `periods.py`,
  `periods.json`, `corrections.py`, `corrections.json` (new), `expected_output.txt`). Both records are re-tracked with a
  background renewed every second, which follows a camera shift at 52.6 to 54.5 s of IMG_8754; bench-v1.0 lost the
  steel sphere there and used 17 of its swings against 20 of the magnets, so the steel period came from earlier, wider
  swings. Both bobs now use the same 20 swings: the period difference is 0.80 +- 0.03 % (bench-v1.0: 0.758 %), and 0.47
  to 0.53 % +- 0.06 % after the corrections (bench-v1.0: 0.41 to 0.48 % +- 0.12 %). The decay correction is now measured
  from the swing amplitudes (bench-v1.0: an estimate, 0.10 +- 0.10 %). The length from the pivot is 1.98 m (bench-v1.0:
  1.96 m). `corrections.py` prints the conditional mass-ratio equivalent, 0.989 to 0.991 (bench-v1.0: 0.9905 to 0.9919),
  and the effective-length equivalent, 18.5 to 21.2 mm.
- **Floating magnet** (`track_waterline.py`, `tracks/IMG_8765_plastic_track.csv`, `tracks/IMG_8766_plastic_track.csv`
  (new); `fit_drift.py`, `expected_output.txt`; `SHA256SUMS_videos.txt`). Two runs on a plastic support are added with
  their tracker, tracks, fits and video hashes; the six earlier fits are unchanged. The basement runs on a wooden stand
  (IMG_8750, IMG_8751) are withdrawn, because the stand held steel screws: their tracks, their tracker
  `track_template.py` and their video hashes are removed.
- **Spinning disc** (`fit_law.py`, `fit_law.json` (new); `expected_output.txt`). The principal fit uses the four
  readings clear of the stop; the seven-reading fit of bench-v1.0 is kept as a sensitivity fit. The readings are
  unchanged; their fits are appended to the expected output.
- **Figures** (`make_figures.py`). Fig. 10 (a) shows the extra time the magnets take to reach each turning point and
  reads `corrections.json`; Fig. 12 (b) shows the two runs on the plastic support instead of the withdrawn basement run;
  Fig. 13 draws the principal and sensitivity fits from `fit_law.json`, with the uncertain
  readings open.
- **Galaxies** (`galaxy_fit.py`, `galaxies.csv`, `expected_output.txt`; `summary.py`, equation number only). G is
  aligned with the paper (6.67408e-11 instead of 6.671e-11). Fitted lambda_g changes in 4 of the 497 curves and m in
  3; the mean standard deviation of the residuals stays 8031 m/s, the mean rms residual goes from 10429 to 10430 m/s,
  the Newtonian-limit mean from 67782 to 67795 m/s; the arctangent baseline (6255 m/s) is unchanged.
  `milky_way.py`, `milky_way_sofue2020.csv` and `milky_way.json` are new: they fit the Milky Way value of the paper,
  which bench-v1.0 did not contain; their output is appended to `expected_output.txt`.
- **Free fall** (`ff_refit.py`, `expected_output.txt`). The script prints its output paths relative to this folder; no
  number changes.

## Licence

See `LICENSE.md`: code under GPL-3.0, data under CC BY-NC 4.0.

## Citation

P. Bouchard, *Finite Theory bench experiments: records and analysis scripts*, release bench-v1.1 (2026),
https://github.com/philippeb8/finite-theory/blob/master/bench-v1.1.zip
