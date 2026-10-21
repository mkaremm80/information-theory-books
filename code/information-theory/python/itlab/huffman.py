"""huffman.py -- Huffman coding (Lab 3).

Conventions (needed to get the same codewords every time)
  * the priority queue holds (probability, order, node); 'order' is a counter
    that breaks ties: older nodes come first, leaves are entered in sorted
    symbol order;
  * at each step the two least probable nodes are removed; the FIRST removed
    becomes the child with bit '0', the second the child with bit '1'.
"""
import heapq
from itertools import count
from math import log2


def huffman(probs):
    """probs: dict symbol -> probability (or count). Returns dict symbol -> code."""
    if len(probs) == 1:
        return {next(iter(probs)): "0"}
    tick = count()
    heap = [(p, next(tick), s) for s, p in sorted(probs.items())]
    heapq.heapify(heap)
    merges = []                                    # for the step tables
    while len(heap) > 1:
        p0, _, n0 = heapq.heappop(heap)            # least probable  -> '0'
        p1, _, n1 = heapq.heappop(heap)            # second least    -> '1'
        merges.append((n0, p0, n1, p1, p0 + p1))
        heapq.heappush(heap, (p0 + p1, next(tick), (n0, n1)))
    huffman.merges = merges
    codes = {}

    def walk(node, prefix):
        if isinstance(node, tuple):                # internal node
            walk(node[0], prefix + "0")
            walk(node[1], prefix + "1")
        else:
            codes[node] = prefix
    walk(heap[0][2], "")
    return codes


def encode(message, codes):
    return "".join(codes[s] for s in message)


def decode(bits, codes):
    inverse = {c: s for s, c in codes.items()}
    out, cur = [], ""
    for b in bits:
        cur += b
        if cur in inverse:
            out.append(inverse[cur])
            cur = ""
    if cur:
        raise ValueError("bitstream ended inside a codeword")
    return "".join(out)


def average_length(probs, codes):
    total = sum(probs.values())
    return sum(p * len(codes[s]) for s, p in probs.items()) / total


def leaves(node):
    return [node] if not isinstance(node, tuple) else leaves(node[0]) + leaves(node[1])


if __name__ == "__main__":
    probs = {"A": 0.30, "B": 0.25, "C": 0.15, "D": 0.12, "E": 0.10, "F": 0.08}
    codes = huffman(probs)
    for step, (n0, p0, n1, p1, p) in enumerate(huffman.merges, 1):
        print(f"step {step}: merge {''.join(leaves(n0))}({p0:.2f}) + "
              f"{''.join(leaves(n1))}({p1:.2f}) = {p:.2f}")
    for s in sorted(codes):
        print(s, codes[s])
    H = -sum(p * log2(p) for p in probs.values())
    L = average_length(probs, codes)
    print("H = %.4f  L = %.4f  efficiency = %.2f%%" % (H, L, 100 * H / L))
    msg = "ABACADEF"
    bits = encode(msg, codes)
    print("encode", msg, "->", bits, "->", decode(bits, codes))
