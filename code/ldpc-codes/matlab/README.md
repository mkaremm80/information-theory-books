# LDPC Codes: From Zero to Expert (Part Two) - MATLAB code (Appendix C)

Put all `.m` files, `5G_bg1.csv`, `5G_bg2.csv` and `pyref_ldpc.mat` in one folder and run:

* `test_ldpc_matlab` - compares every function with the Python reference results (10 tests)
* `demo_ldpc`        - worked examples of the book
* `sim_ldpc_matlab`  - FER of the 5G NR code K = 1024, R = 1/2
* `final_project_lab`- transport block with HARQ (rv0 + rv2)

MATLAB R2016b or newer, no toolbox. Tested with GNU Octave 8.4.
Indices are 1-based; frames are rows (B x n).
