"""Dataset types and interval helpers."""

import numpy as np
import pytest

from dadl.data import ChannelId, Dataset, TagSet, tags_sha256
from dadl.stats import Laplace, beta_posterior, clopper_pearson, gamma_posterior, wilson


def test_channel_id_validates_and_orders():
    assert ChannelId("a", 1) < ChannelId("a", 2) < ChannelId("b", 0)
    with pytest.raises(ValueError):
        ChannelId("", 0)
    with pytest.raises(ValueError):
        ChannelId("a", -1)
    with pytest.raises(TypeError):
        ChannelId("a", 1.5)


def test_tagset_is_validated_and_read_only():
    ch = ChannelId("dev", 0)
    ts = TagSet({ch: np.array([1, 5, 9])}, clock="sim")
    with pytest.raises(ValueError):
        ts[ch][0] = 7
    with pytest.raises(TypeError):
        ts.tags[ch] = np.array([1])
    assert ts.duration(ch) == 8
    assert TagSet({ch: np.array([1, 5])}, duration_ps=100).duration() == 100
    with pytest.raises(ValueError, match="sorted"):
        TagSet({ch: np.array([5, 1])})
    with pytest.raises(TypeError):
        TagSet({ch: np.array([0.5, 1.0])})
    with pytest.raises(TypeError):
        TagSet({"A+": np.array([1])})


def test_dataset_backend():
    with pytest.raises(ValueError):
        Dataset({}, backend="lab")


def test_tag_hash_depends_on_values_not_dtype():
    assert tags_sha256(np.array([1, 2], dtype=np.int32)) == tags_sha256(np.array([1, 2], dtype=np.int64))
    assert tags_sha256(np.array([1, 2])) != tags_sha256(np.array([1, 3]))


def test_wilson_and_clopper_pearson_known_values():
    w = wilson(10, 100)
    assert (round(w.low, 4), round(w.high, 4)) == (0.0552, 0.1744)
    cp = clopper_pearson(0, 10)
    assert cp.low == 0.0 and round(cp.high, 4) == 0.3085
    with pytest.raises(ValueError):
        wilson(5, 3)


def test_conjugate_posteriors():
    rng = np.random.default_rng(0)
    iv, draws = gamma_posterior(400, 2.0, rng, 4000)
    assert iv.value == 200 and iv.contains(200)
    assert np.quantile(draws, 0.5) == pytest.approx(200, rel=0.02)
    iv, draws = beta_posterior(30, 100, rng, 4000)
    assert iv.contains(0.3) and np.all((draws > 0) & (draws < 1))


def test_laplace_refuses_an_unidentified_fit():
    with pytest.raises(ValueError, match="not identified"):
        Laplace.from_hessian(np.zeros(2), np.array([[1.0, 0.0], [0.0, -1.0]]))
