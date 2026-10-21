"""Write ../ldpc_listing.tex: ldpc.py printed in segments (one box per block)."""
import re
lines = open("ldpc.py").read().splitlines()
starts = [1] + [i for i, l in enumerate(lines, 1) if re.match(r"# ==== \[\w+\]", l) and i > 1]
starts = sorted(set(starts))
segs = []
for a, b in zip(starts, starts[1:] + [len(lines) + 1]):
    segs.append((a, b - 1))
out = []
for a, b in segs:
    out.append(r"\pyfile[listing options={style=py,firstline=%d,lastline=%d,firstnumber=%d},"
               r"title={\ttfamily\small ldpc.py \textbullet\ lines %d--%d},fonttitle=\color{white},colbacktitle=navy3]{code/ldpc.py}" % (a, b, a, a, b))
open("../ldpc_listing.tex", "w").write("\n".join(out) + "\n")
print(segs)
