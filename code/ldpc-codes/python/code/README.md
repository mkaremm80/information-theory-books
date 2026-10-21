# LDPC Codes: From Zero to Expert (Part Two) - Python code

Library: `ldpc.py` (NumPy only). Scripts:

| file | purpose |
|---|---|
| `test_ldpc.py` | 10 unit tests - run this first |
| `examples.py` | every number of the worked examples (Chapters 2-11) |
| `demo.py` | demonstration of Chapter 12 |
| `solved.py` | numbers of the solved problems (Chapter 14) |
| `run_sims.py A|B|C|D` | Monte-Carlo runs of Chapter 13 |
| `compare_polar.py` | LDPC versus polar (uses `polar.py` from Part One) |
| `design_bec.py` | linear-programming design of irregular codes (Chapter 10, needs SciPy) |
| `project.py` | final project: 5G NR transport block with HARQ (Chapter 16) |
| `make_figs.py`, `plot_sims.py` | figures (Matplotlib) |
| `make_pyref.py` | reference data for the MATLAB tests |

The folder `../data` must contain `5G_bg1.csv` and `5G_bg2.csv` (base-graph tables of TS 38.212, taken from the
open-source Sionna library, Apache-2.0 licence; the licence file is included).
