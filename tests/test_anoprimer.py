import numpy as np
import pandas as pd
import pytest

import AnoPrimer

ace1_seq = "GCGGCGGCTTCTACTCCGG"
kdr_seq = "AGTGATAGGAAATTTAGTCGT"


@pytest.mark.parametrize(
    "sequence",
    [ace1_seq, kdr_seq],
)
def test_check_my_oligo(sequence):
    AnoPrimer.check_my_oligo(sequence=sequence, sample_sets="AG1000G-GH")
