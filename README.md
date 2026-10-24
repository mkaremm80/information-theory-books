# Information Theory Books

Books on information theory and channel coding, with companion MATLAB and Python code.

By **Mahmoud Abdelaziz** — Associate Professor, Avionics Department, MTC; IEEE Senior Member.

Live showcase: https://mkaremm80.github.io/information-theory-books/


## ▶ Watch: 5G NR LDPC in 60 seconds

**How the encoder builds a codeword — and how the layered decoder repairs it.**

<video controls src="videos/ldpc-5g/ldpc5g_encoder_decoder_v1.mp4" poster="videos/ldpc-5g/poster.jpg"></video>

[▶ Watch the video (60 s, MP4)](videos/ldpc-5g/ldpc5g_encoder_decoder_v1.mp4) · [📄 Narration, timeline and numbers](videos/ldpc-5g/ldpc5g_encoder_decoder_v1_narration.md) · [📘 LDPC book (Part 2)](books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf)

One real 5G code block: **K = 1024 data bits, E = 2048 sent bits, rate 1/2, base graph BG2, lifting size Z = 104, BPSK at Eb/N0 = 1.5 dB, layered normalized min-sum (α = 0.75).**

| | At the receiver | After 10 iterations |
|---|---|---|
| Unsatisfied checks (of 1248) | 620 | **0** |
| Wrong bits (of 2288) | 373 | **0** |

Over many frames this code has FER ≈ 4 × 10⁻³ at 1.5 dB, with 9.7 iterations on average. From *LDPC Codes: From Zero to Expert (Part 2)*, Chapters 6, 9 and 11.


## Books

| Book | PDF | Code |
|---|---|---|
| Information Theory: Compression | [books/information-theory/Information_Theory_Compression_v1.pdf](books/information-theory/Information_Theory_Compression_v1.pdf) | [code/information-theory/](code/information-theory/) |
| LDPC Codes: From Zero to Expert (Part 2) | [books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf](books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf) | [code/ldpc-codes/](code/ldpc-codes/) |
| Polar Codes: From Zero to Expert | [books/polar-codes/Polar_Codes_From_Zero_to_Expert_v2.pdf](books/polar-codes/Polar_Codes_From_Zero_to_Expert_v2.pdf) | [code/polar-codes/](code/polar-codes/) |

For Polar Codes, version 2 (October 2026, 156 pages) is now one volume with block diagrams and Appendix D (MATLAB) included inside it, so there is no separate MATLAB appendix file. The first version (v1) and a draft of chapters 1–9 are also in the book folder. The Information Theory book folder still includes its MATLAB appendix.

## Code

- **MATLAB** (`.m` files): entropy, Huffman, arithmetic and LZW coding, Hamming/BCH codes, Markov sources, LDPC codes (belief propagation, 5G base graphs), and polar codes (SC/list decoding, CRC-aided decoding).
- **Python** (`.py` files): the same topics as a Python package, with test scripts.

Run the test scripts (`test_*.m`, `test_*.py`) or the demo scripts (`demo_*.m`, `demo.py`) to verify the examples.

## Note

The book PDFs live under `books/`, and the code folders also include the MATLAB (`.mat`) and NumPy (`.npy`) data files the scripts use.

## License

- Code in `code/` is released under the MIT License (see [LICENSE](LICENSE)).
- Books in `books/` are licensed under Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) (see [CC-BY-NC-4.0.txt](CC-BY-NC-4.0.txt)): free to share and adapt with attribution, for non-commercial purposes only.
