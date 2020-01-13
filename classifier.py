import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from bioinfenhancers.transfomer import Transfomer

from Bio import SeqIO

positive_records = list(SeqIO.parse("data/vista1500", "fasta"))
negative_records = list(SeqIO.parse("data/randoms1500", "fasta"))


positive_data = [(str(record.seq).upper(), 1) for record in positive_records]
negative_data = [(str(record.seq).upper(), 0) for record in negative_records]


data = positive_data + negative_data
data[:5]

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

