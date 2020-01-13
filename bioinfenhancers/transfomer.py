from Bio.Seq import Seq
import itertools
import mmh3
from typing import List, Tuple, Dict


def chunkstring(string, length, step):
    return (string[0 + i : length + i] for i in range(0, len(string), step))


def hash_kmer(kmer):
    rc_kmer = str(Seq(kmer).reverse_complement())

    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer

    return canonical_kmer


class Transfomer:
    def __init__(self, k: int = 4):
        alphabet = "ACTG"
        self._combinations = [
            "".join(output) for output in itertools.product(alphabet, repeat=k)
        ]

    def get_counter(self) -> Dict[str, int]:
        return {hash_kmer(x): 0 for x in self._combinations}

    def get_feature_vector(self, frame: str) -> Tuple[bool, List[float]]:
        counter = {hash_kmer(x): 0 for x in self._combinations}
        valid = True
        for quad in chunkstring(frame, 4, 1):
            if "N" in quad:
                valid = False
                break
            if len(quad) == 4:
                counter[hash_kmer(quad)] += 1

        counts = [item for key, item in sorted(counter.items())]
        freq = [x / 1500 for x in counts]
        return valid, freq
