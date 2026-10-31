# Finite Theory bench experiments: records and analysis scripts

Release **bench-v1.0** (29 October 2026). These are the records and scripts behind Section 5 (Laboratory experiments) of
*Finite Theory: Solutions for Physics Anomalies*, P. Bouchard, Services Informatiques Fornux. Every number in Tables 3 to 6
and Figures 9, 10, 12 and 13 of the paper is reproduced by the commands below. Each folder holds an `expected_output.txt`
with the output of this release.

## Contents

| Folder | Experiment | Records | Scripts | Paper |
|---|---|---|---|---|
| `freefall/` | 40 paired drops, NdFeB and steel spheres | `data/day7..day11/FreeFall-7pm-{Up,Down}.py` (tracked positions, 240 fps, in their `DATA` blocks) | `ff_refit.py` | Sec. 5.1, Table 3, Fig. 9 |
| `pendulum/` | Magnets and steel sphere on the same nylon twine | `tracks/IMG_8754_steel.csv`, `tracks/IMG_8755_magnets.csv` | `track_pendulum.py`, `periods.py`, `corrections.py` | Sec. 5.2, Table 4, Fig. 10 |
| `floating_magnet/` | Float in the afternoon, at night and in the basement; phone magnetometer | `tracks/*.csv`, `magnetometer/*.csv` | `track_night.py`, `track_template.py`, `fit_drift.py`, `magnetometer/holds.py` | Sec. 5.3, Table 5, Figs. 11 and 12 |
| `spinning_disc/` | Disc on a beam balance | readings in `readings.py` | `readings.py` | Sec. 5.4, Table 6, Fig. 13 |
| `figures/` | | | `make_figures.py` (writes `figures/out/`) | Figs. 9, 10, 12, 13 |
| `galaxies/` | Rotation-curve survey: 497 curves of SPARC and the Sofue atlas | read from github.com/philippeb8/grc at commit 9f3c865 | `galaxy_fit.py`, `baseline_arctan.py`, `summary.py`; fitted parameters in `galaxies.csv` | Sec. 4.3.1 |
| `archive/` | Interferometer experiment proposed earlier for the ISS, with its open steps | | `ISS_proposal_archived.tex` and `.pdf` | removed from the paper |

The gallium-vortex run (Sec. 5.5) is a single reading, 2.6 +- 3.7 g, and has no file.

## Run

Python 3.10 or later with numpy, scipy and matplotlib (`pip install -r requirements.txt`). From this folder:

```
python3 freefall/ff_refit.py freefall/data       # accelerations by half of the drop, Table 3; writes freefall/lower_half.json
python3 pendulum/periods.py                       # periods of both runs; writes pendulum/periods.json
python3 pendulum/corrections.py                   # Table 4
python3 floating_magnet/fit_drift.py              # Table 5
python3 floating_magnet/magnetometer/holds.py     # field at each 10 cm hold
python3 spinning_disc/readings.py                 # Table 6 readings and calibration
python3 figures/make_figures.py                   # Figs. 9, 10, 12, 13 (needs the two json files above)
```

The galaxy survey needs the data of the GRC repository and takes a few minutes:

```
git clone https://github.com/philippeb8/grc && git -C grc checkout 9f3c8658fd7be9aa43964caf29e8ed39a0d4f575
cd galaxies && python3 galaxy_fit.py ../grc && python3 baseline_arctan.py ../grc && python3 summary.py
```

Re-tracking from the videos needs ffmpeg: `track_pendulum.py <video> <out.csv>`, `track_night.py <video> <name>`,
`track_template.py <video> <out.csv> <x0> <x1>`. The videos are available from the author on request;
`SHA256SUMS_videos.txt` identifies the files the tracks were made from.

## Notes on the records

- **Free fall.** `ff_refit.py` reads the `DATA` blocks of the ten session scripts without executing them. The per-session
  summary that those scripts print fits the median of each block with the release height fixed at zero; it is biased and
  is not used. Each drop is calibrated on its own steel sphere; the drop with a failed steel track (day 8, D) is excluded.
  `per_drop.csv` lists, for every drop, the calibration scale and curvature, the residuals of the calibration and of the
  two-half fit, the mid-height time and the fitted accelerations.
- **Pendulum.** The correction for the slower decay of the magnets' swing (0.10 +- 0.10 %) is an estimate, and the stretch
  of the nylon twine is estimated from the modulus of nylon (0.5 to 2 GPa); both are marked in `corrections.py`.
- **Floating magnet.** The afternoon tracks (IMG_8733, IMG_8735) were made with an earlier version of the
  background-subtraction tracker, which was not kept; the tracks are the record. All six drifts are fitted by the same
  script. The basement run IMG_8750 tracks poorly and is described from its frames in the paper.
- **Spinning disc.** The pointer positions were read from the videos; `readings.py` holds them with the coin calibration.
- **Galaxies.** `galaxy_fit.py` follows `src/main.cpp` of the GRC repository: lambda_g on a grid from 1e18 to 1e22 kg/m, m between
  1/15 and 15 times a photometric mass scale, omega (rad/s, the angular velocity of the comoving framework of the galaxy)
  set so the curve passes through the last point, objective = standard
  deviation of the residuals about their mean (unweighted). `baseline_arctan.py` fits an empirical arctangent curve with the
  same linear term, objective and parameter count. Duplicate objects are kept, as in the original survey, and listed.

## Licence

See `LICENSE.md`: code under GPL-3.0, data under CC BY-NC 4.0.

## Citation

P. Bouchard, *Finite Theory bench experiments: records and analysis scripts*, release bench-v1.0 (2026),
https://github.com/philippeb8/finite-theory/tree/main/bench
