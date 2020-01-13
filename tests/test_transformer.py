from bioinfenhancers.transfomer import hash_kmer, Transfomer


def test_same_hash():
    assert hash_kmer("ATGGC") == hash_kmer("GCCAT")


def test_combinations_count():
    transformer = Transfomer(k=4)
    assert len(transformer.get_counter().keys()) == 136


def test_invalid_sequence():
    transformer = Transfomer(k=4)
    valid, _ = transformer.get_feature_vector("AAAANNNN")
    assert valid == False
