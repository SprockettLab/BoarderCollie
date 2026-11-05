from boardercollie.repeats import reverse_complement, unit_identity


def test_identical_units():
    assert unit_identity("ACGTTGCA", "ACGTTGCA", 0.1) == 1.0


def test_rotation_is_identical():
    assert unit_identity("ACGTTGCA", "TGCAACGT", 0.1) == 1.0


def test_reverse_complement_is_identical():
    u = "AACCGGTTAGC"
    assert unit_identity(u, reverse_complement(u), 0.1) == 1.0


def test_one_mismatch_in_sixteen():
    a = "ACGTTGCAAGGCTTAC"
    b = "ACGTTGCAAGGCTTAA"
    assert abs(unit_identity(a, b, 0.1) - 15 / 16) < 1e-9


def test_length_outside_tolerance_is_zero():
    assert unit_identity("ACGTACGTAC", "ACGTACGTACGTACGT", 0.1) == 0.0


def test_empty_unit_is_zero():
    assert unit_identity("", "ACGT", 0.1) == 0.0
