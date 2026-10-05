# ftsim — N point charges under Finite Theory's time-dilation factor

## Build

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/ftsim
```

Requires Qt6 Widgets (a Qt5 fallback is commented in `CMakeLists.txt`).
`core_test.cpp` is a standalone check of the acceleration kernel and needs no Qt:

```sh
g++ -O2 -std=c++17 core_test.cpp -o core_test && ./core_test
```

## Physics

`f(r) = eta*r / (eta*r + S)`, with `S = |q_source|`.

| route | acceleration | behaviour |
|---|---|---|
| Potential | `-C·eta²(S − eta·r)/(eta·r + S)³` | sign change at `r = S/eta` |
| Centripetal | `C·eta²/(eta·r + S)²` | sign-definite, no crossover |
| Coulomb | `C/r²` | reference |

`C = K·q_i·q_j/m_i`, and `K = 1` in the unit system below.

Both FT routes are **finite at r = 0**, so the integrator needs no softening
length — the `f` factor regularises the Coulomb singularity by itself.

## Units

| quantity | unit | value |
|---|---|---|
| length | a0 | 5.29177210903e-11 m |
| speed | v_Bohr | 2.18769126e6 m/s |
| time | a0/v_Bohr | 2.4189e-17 s |
| charge | e | 1.602176634e-19 C |
| mass | m_e | 9.1093837015e-31 kg |

`ETA = eta·a0/e`. **ETA = 8.7705e5** corresponds to `eta = c·eps0`; the
crossover is then at 1.14e-6 a0 and the run is indistinguishable from Coulomb
(ratio to Coulomb at 1 a0 = 0.999995439). Set ETA to roughly 0.1–10 to make
the FT structure visible on screen.

## License

Copyright (C) 2026 Phil Bouchard <phil@fornux.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See `COPYING` for the full text.
