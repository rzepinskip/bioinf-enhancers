from bioinfenhancers.transfomer import hash_kmer, Transfomer


def test_same_hash():
    assert hash_kmer("ATGGC") == hash_kmer("GCCAT")


def test_combinations_count():
    transformer = Transfomer(k=4)
    assert len(transformer.get_counter().keys()) == 136
