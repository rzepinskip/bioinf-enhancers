import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from bioinfenhancers.transfomer import Transfomer
from Bio import SeqIO

positive_records = SeqIO.parse("data/vista1500", "fasta")
negative_records = SeqIO.parse("data/randoms1500", "fasta")
positive_data = [(str(record.seq).upper(), 1) for record in positive_records]
negative_data = [(str(record.seq).upper(), 0) for record in negative_records]
data = positive_data + negative_data

transfomer = Transfomer(k=4)
dataset = np.array(
    [transfomer.get_feature_vector(frame)[1] + [label] for frame, label in data]
)

X = dataset[:, 0:136]
y = dataset[:, 136]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=44
)

model = RandomForestClassifier()
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

predictions = [round(value) for value in y_pred]

accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))

chr_data = pd.read_csv("frequencies.csv")

chr_X = chr_data.iloc[:, 2:].values

chr_y_pred = [true_prob for _, true_prob in model.predict_proba(chr_X)]

chr_data["prob"] = chr_y_pred
avg_valid_prob = chr_data[chr_data.valid == 1].prob.mean()
chr_data.loc[chr_data.valid == 0, "prob"] = avg_valid_prob

# chr_data.sample(10)[["index", "valid", "prob"]]

with open("data/chr21.wig", "w") as f:
    f.write("fixedStep chrom=chr21 start=0 step=750 span=1500\n")
    f.write("\n".join([str(x) for x in chr_data.prob.tolist()]))
