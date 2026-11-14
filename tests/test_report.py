import math

from orf_profiler.report import median_from_counts


def test_median_odd_and_even():
    assert median_from_counts({10: 1, 20: 1, 30: 1}) == 20
    assert median_from_counts({10: 2, 20: 2}) == 15
    assert median_from_counts({10: 5}) == 10


def test_median_weighted():
    assert median_from_counts({10: 1, 50: 9}) == 50


def test_median_of_nothing_is_nan():
    assert math.isnan(median_from_counts({}))
