"""Write ../linemap.tex: first/last line of every [block] of ldpc.py."""
import re
lines = open("ldpc.py").read().splitlines()
out = []
for i, l in enumerate(lines, 1):
    m = re.match(r"# ==== \[(/?)(\w+)\]", l)
    if m and not m.group(1):
        out.append(rf"\expandafter\def\csname lm@{m.group(2)}@a\endcsname{{{i + 1}}}")
    elif m:
        out.append(rf"\expandafter\def\csname lm@{m.group(2)}@b\endcsname{{{i - 1}}}")
open("../linemap.tex", "w").write("\n".join(out) + "\n")
print("\n".join(out))
