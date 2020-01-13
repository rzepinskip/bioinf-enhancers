import pandas as pd
import gzip
import logging

from Bio import SeqIO
from bioinfenhancers.transfomer import chunkstring, Transfomer

logging.basicConfig(level=logging.DEBUG)

with gzip.open("data/chr21.fa.gz", "rt") as handle:
    record = next(SeqIO.parse(handle, "fasta"))
frames = list(chunkstring(str(record.seq).upper(), 1500, 750))
transformer = Transfomer(k=4)

output = []
for index, frame in enumerate(frames):
    if index % 1000 == 0:
        logging.debug(f"{index+1}/{len(frames)}")
    valid, freq = transformer.get_feature_vector(frame)
    output.append([index, int(valid)] + freq)

df = pd.DataFrame(
    output,
    columns=["index", "valid"]
    + [key for key, item in sorted(transformer.get_counter().items())],
)

df.to_csv("frequencies.csv", index=False)
