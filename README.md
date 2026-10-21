# Information Theory Books

Books on information theory and channel coding, with companion MATLAB and Python code.

By **Mahmoud Abdelaziz** — Associate Professor, Avionics Department, MTC; IEEE Senior Member.

## Books

| Book | PDF | Code |
|---|---|---|
| Information Theory: Compression | [books/information-theory/Information_Theory_Compression_v1.pdf](books/information-theory/Information_Theory_Compression_v1.pdf) | [code/information-theory/](code/information-theory/) |
| LDPC Codes: From Zero to Expert (Part 2) | [books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf](books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf) | [code/ldpc-codes/](code/ldpc-codes/) |
| Polar Codes: From Zero to Expert | [books/polar-codes/Polar_Codes_From_Zero_to_Expert_v2.pdf](books/polar-codes/Polar_Codes_From_Zero_to_Expert_v2.pdf) | [code/polar-codes/](code/polar-codes/) |

The book folders also include appendices (MATLAB appendices) and, for Polar Codes, the first version (v1) and a draft of chapters 1–9.

## Code

- **MATLAB** (`.m` files): entropy, Huffman, arithmetic and LZW coding, Hamming/BCH codes, Markov sources, LDPC codes (belief propagation, 5G base graphs), and polar codes (SC/list decoding, CRC-aided decoding).
- **Python** (`.py` files): the same topics as a Python package, with test scripts.

Run the test scripts (`test_*.m`, `test_*.py`) or the demo scripts (`demo_*.m`, `demo.py`) to verify the examples.

## Note

The book PDFs live under `books/`, and the code folders also include the MATLAB (`.mat`) and NumPy (`.npy`) data files the scripts use.
