# Polar codes — polarization, encoder and SC decoder: explainer video (72 s)

▶ [Watch the video](polar_codes_explainer_v1.mp4)

**In 72 seconds:** how channels polarize, how the encoder works, and how the SC decoder uses it.
Start with BEC(0.5): at N = 1024, 34% of bit-channels are almost perfect and 34% almost useless, mean capacity still 0.5. Then a small (8,4) code, and a (256,128) result: SC BLER ≈ 1.2×10⁻² at 3 dB, CA-SCL-8 ≈ 6.6×10⁻⁴.

📘 Source book: [Polar Codes: From Zero to Expert](../../books/polar-codes/Polar_Codes_From_Zero_to_Expert_v2.pdf), Chapters 4 (polarization), 5 (encoding), 6 (SC decoding) and 7–8 (list decoding results).

---
Source: *Polar Codes: From Zero to Expert* (v2), Chapters 4 (polarization), 5 (encoding), 6 (SC decoding) and 7–8 (list decoding results).
All numbers come from the book's Python library (`polar.py`): exact BEC(0.5) polarization, the (8,4) code of Worked Examples 5.3 and 6.3, and the (256,128) Monte-Carlo results of Fig. 7.3.

## Timeline

| time | caption | what to say (simple English) |
|---|---|---|
| 0–4 s | 1 · Polar codes are built from one idea… | "Polar codes come from one idea by Arıkan, in 2009. We combine N copies of a channel, then we split them again into N bit-channels. We start with an erasure channel that loses half of the bits. Its capacity is 0.5." |
| 4–9 s | 2 · One polarization step | "One step: two copies become a worse channel, W minus, and a better channel, W plus. Their erasure probabilities are 0.75 and 0.25. Nothing is lost: 0.25 plus 0.75 is still two times 0.5. The information only moves." |
| 9–14 s | 3 · Repeat the step | "Now we repeat the step on every channel. Look at the tree: after four levels, with 16 channels, many values are already close to 0 or to 1. Orange means useless, blue means almost perfect." |
| 14–26 s | 4 · Grow to N = 1024 | "We keep going up to 1024 channels. The histogram shows two peaks, at 0 and at 1. At N = 1024, 34 percent of the channels are almost perfect and 34 percent are almost useless. The average is still exactly 0.5. A polar code sends data only on the good channels." |
| 26–29 s | 5 · Back to a small code | "Now a small example with 8 bits. The 4 best bit-channels, 3, 5, 6 and 7, carry the message 1 0 1 1. The other 4 are frozen to zero, and the receiver knows them." |
| 29–32 s | 6 · Encoder stage 1 | "The encoder has 3 butterfly stages. In stage 1, the upper half is XORed with the lower half." |
| 32–35 s | 7 · Stages 2 and 3 | "Stages 2 and 3 do the same inside smaller blocks. That is only 12 XOR operations." |
| 35–38 s | 8 · Codeword | "The codeword is 1 0 1 0 0 1 0 1. Each stage is one level of the polarization tree." |
| 38–42 s | 9 · Channel | "We send the bits with BPSK and add noise. The receiver computes the LLRs. Position 4 has the wrong sign, so a simple hard decision makes one error." |
| 42–53 s | 10 · SC decoder, left half | "The SC decoder sends the LLRs down the tree. The f function goes left, the g function goes right. At each leaf it decides one bit: frozen bits are zero, information bits follow the sign of the LLR. The decisions go back up as partial sums." |
| 53–64 s | 11 · Right half | "In the right half, the g function uses the bits that are already known. The LLRs become large: 5.1, 3.5 and 12.2. The good bit-channels really are good. This is polarization working inside the decoder." |
| 64–72 s | 12 · Done | "The decoder finds the correct message, 1 0 1 1, and corrects the channel error. On a longer code with 256 bits, SC has a block error rate of about 1 in 100 at 3 dB, and CRC-aided list decoding with 8 paths is about 18 times better." |

## Panel guide

- **Left (scene):** polarization tree for BEC(0.5), levels N = 1 … 16. Upper branch = minus (worse), lower = plus (better), in natural index order. Colour = capacity (orange 0 → blue 1); the current level has gold frames. Then: butterfly encoder of the (8,4) code, stage by stage, with the bit values after each stage, then BPSK and the received LLRs (red = wrong sign). Then: the SC decoding tree with the LLR vector of each node, the partial sums β (green) and the leaf decisions (gold = current step, 14 steps).
- **Top middle:** polarization phase: capacities of all N bit-channels in natural order. Code phase: the 8 capacities with the information set (blue) and the frozen set (grey).
- **Top right:** polarization phase: histogram of the capacities. Encoder phase: the 8 bits after each stage. Decoder phase: channel LLRs α (teal, red = wrong sign) next to the decision LLRs of the bit-channels (blue = information, grey = frozen).
- **Bottom middle:** polarization phase: fraction of almost perfect / almost useless / not yet polarized channels up to n = 20 (book Fig. 4.6, thresholds Z < 10⁻³ and Z > 1 − 10⁻³); the red dots move with n. For ε = 0.5 the blue and orange curves are identical (symmetry). Code phase: BLER of the (256,128) code with SC, SCL (L = 8) and CA-SCL (L = 8, CRC-11) (book Fig. 7.3).
- **Bottom right:** polarization phase: capacities sorted from best to worst; the best half is the information set of a rate-1/2 code (mean capacities of both halves shown). Decoder phase: table of the SC decisions (set, sent bit, decision LLR, decision, ✓).

## Numbers in the video (checked by the script)

| quantity | value |
|---|---|
| encoder stages (u → x) | 00010011 → 00100011 → 10101111 → 10100101 |
| leaf LLRs L of bit-channels 0…7 | −0.6, +0.6, +0.2, −3.4, +0.9, +5.1, −3.5, −12.2 |
| decisions û | 0 0 0 1 0 0 1 1 (correct) |
| N = 1024, BEC(0.5) | 34% Z < 10⁻³, 34% Z > 1 − 10⁻³, mean capacity 0.5 |
| (256,128) BLER at 2 dB | SC 1.3×10⁻¹, SCL-8 4.2×10⁻², CA-SCL-8 4.0×10⁻² |
| (256,128) BLER at 3 dB | SC 1.2×10⁻², SCL-8 3.7×10⁻³, CA-SCL-8 6.6×10⁻⁴ |

## Settings

BEC(0.5), n swept 1 → 10 (N = 2 … 1024) with an eased schedule; the (8,4) code with F = {0,1,2,4}, m = 1011; received LLRs α = (−1.5, 2.0, −0.8, 1.2, −0.6, −2.5, 1.8, −3.0) exactly as in Worked Example 6.3; min-sum f. Timings: intro 4 s, polarization 22 s, encoder 12 s, channel 4 s, decoder 22 s, final hold 8 s; 1920 × 1080, 30 fps.

## Easy edits (top of `polar_video.py`)

- `EPS`: try 0.3 (good channels win) or 0.7; the curves in the bottom-middle panel then separate.
- `ALPHA`: put another noisy LLR vector, for example a stronger error, to show an SC failure.
- `FROZEN`, `U`: another (8,K) code; keep `U[FROZEN] = 0`.
- `T_POL`, `T_DEC`, …: timings in seconds.

Run `python3 polar_video.py --preview` (six PNG frames) before `python3 polar_video.py --out video.mp4`. The script needs `polar.py` and `data/sim_A.json` from the book's code archive (set `BOOK` at the top).