# 5G NR LDPC encoder and decoder — explainer video (60 s)

▶ [Watch the video](ldpc5g_encoder_decoder_v1.mp4)

**In 60 seconds:** how the encoder builds a 5G codeword, and how the layered decoder repairs it after the channel.
One block: K = 1024, E = 2048, rate 1/2, BG2, Z = 104, BPSK at 1.5 dB. Receiver starts with 620 unsatisfied checks and 373 wrong bits — after 10 iterations: 0 and 0.

📘 Source book: [LDPC Codes: From Zero to Expert (Part 2)](../../books/ldpc-codes/LDPC_Codes_From_Zero_to_Expert_Part2_v1.pdf), Chapters 6 (encoding), 9 (layered decoding) and 11 (5G NR).

---
Source: *LDPC Codes: From Zero to Expert (Part Two)*, Chapters 6 (encoding), 9 (layered decoding) and 11 (5G NR).
All numbers come from the book's Python library (`ldpc.py`): the real TS 38.212 base graph BG2 and one real noisy frame.

## Timeline

| time | caption | what to say (simple English) |
|---|---|---|
| 0–5 s | 1 · One 5G NR code block … | "Here is one 5G code block: 1024 data bits at rate one half. The rule of TS 38.212 picks base graph 2 with lifting size 104. Every coloured square is a 104 by 104 shifted identity, and the number is its shift. We only need 12 of the 42 block rows, because we send 2048 bits." |
| 5–10 s | 2 · Encoder step 1: λ = A·s | "The encoder works in three steps. Step one: multiply the data bits by the sparse part A. Only the four core rows are used. Sixteen filler bits, known zeros, fill the block up to 1040 bits." |
| 10–16 s | 3 · Encoder step 2: core parity | "Step two: the core. Add the four rows of lambda: the dual-diagonal blocks cancel, and we get the first parity block p-zero. Then back-substitution gives the other three blocks: 416 parity bits." |
| 16–21 s | 4 · Encoder step 3: extension parity | "Step three: each extension row has its own identity block, so each new parity bit is just the XOR of the other bits in its check. Now every one of the 1248 checks is satisfied." |
| 21–25 s | 5 · Channel | "The first 208 bits are punctured: they are never sent. We send 2048 bits with BPSK at 1.5 dB. The receiver starts with 373 wrong hard decisions and 620 broken checks." |
| 25–30 s | 6 · Layered NMS decoder | "The decoder takes one block row at a time — 104 checks in parallel — and updates the bit beliefs at once. Watch the gold band move down the base graph." |
| 30–36 s | 7 · Iterations 1–4 | "The punctured bits start with no information, but their many checks recover them. The red marks disappear, and the histogram of beliefs moves to the right." |
| 36–50 s | 8 · Iterations 5–10 | "Now the errors fall fast. After every iteration the decoder checks the syndrome. When it is zero, it stops." |
| 50–60 s | 9 · Done after 10 iterations | "This frame needed 10 iterations. Over many frames, the book shows this code has a frame error rate of about 4 in 1000 at 1.5 dB, with 9.7 iterations on average." |

## Panel guide

- **Left (scene):** the 12 × 22 part of BG2 that is used (set index 6, shifts P = V mod 104). Grey = punctured columns, blue = systematic, orange = core parity, light orange = extension parity. Gold frame = the encoder step; gold band = the decoder's current layer. The red bars on the right count unsatisfied checks per block row (out of 104).
- **Top middle:** the 2288 code bits, one row per base column. In the encoder phase: dark = 1, grey = filler bits. In the decoder phase: red = wrong hard decision.
- **Top right:** encoder phase: where the bits go (208 punctured + 816 + 416 + 816 sent = 2048). Decoder phase: histogram of the a-posteriori LLRs, multiplied by (1 − 2c) so that positive means correct; the grey bar is LLR = 0 (punctured bits at the start).
- **Bottom middle:** unsatisfied checks, wrong bits and wrong data bits after every layer (12 layers = 1 iteration).
- **Bottom right:** frame error rate of the same code over many frames (book, Fig. 13.3); the red dot is the channel of this frame.

## Numbers in the video (from the run)

| after iteration | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| unsatisfied checks (of 1248) | 620 | 418 | 361 | 298 | 215 | 162 | 67 | 21 | 11 | 1 | 0 |
| wrong bits (of 2288) | 373 | 336 | 297 | 275 | 198 | 152 | 79 | 24 | 14 | 1 | 0 |

The 373 wrong bits at the start include about half of the 208 punctured bits, whose LLR is 0.

## Settings

K = 1024, E = 2048, BG2, Z = 104, set index 6, Eb/N0 = 1.5 dB, layered NMS with α = 0.75, max 50 iterations, seed 1. Seeds 0–39 were scanned; seed 1 converges in 10 iterations, close to the book average of 9.7 at this Eb/N0. (Some seeds converge in 5–7 iterations; seed 31 fails after 50 iterations. Such failures are rare: the book measures FER ≈ 0.4% at this Eb/N0.)

## Easy edits (top of `ldpc5g_video.py`)

- `EBN0`: try 1.25 dB (slower convergence) or 2.0 dB (very fast).
- `SEED`: 31 shows a decoding failure; 21 converges in 5 iterations.
- `K`, `E`: any code block, e.g. `K, E = 1024, 1152` for rate 8/9 (BG1).
- `T_ENC`, `T_DEC`, `T_HOLD`: timings in seconds.

Run `python3 ldpc5g_video.py --preview` (3–4 PNG frames) before `python3 ldpc5g_video.py --out video.mp4`. The script needs `ldpc.py` and the CSV tables from the book's code archive (set `BOOK` at the top).