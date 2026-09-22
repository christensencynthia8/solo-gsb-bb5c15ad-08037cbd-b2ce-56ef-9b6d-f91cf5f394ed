"""部分酶切概率分布接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 线性序列（长度 14）：AAAA GAATTC TTTT，EcoRI 单切口 (top=5, bottom=9)
SEQ = "AAAAGAATTCTTTT"
ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
ECORV = {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}  # SEQ 上无位点

# 线性序列（长度 14）：AAAA GAA A GAC AAA
#   Ava (GAA, 1/1)：平端切口 (5, 5)；Bvo (GAC, 1/1)：平端切口 (9, 9)
SEQ2 = "AAAAGAAAGACAAA"
AVA = {"name": "Ava", "pattern": "GAA", "cut_top": 1, "cut_bottom": 1}
BVO = {"name": "Bvo", "pattern": "GAC", "cut_top": 1, "cut_bottom": 1}

# 环形序列（长度 16）：EcoRI 切口 (1, 5)，EcoRV 切口 (13, 13)
CIRC = "GAATTCGGGGGATATC"


def _create(client, sequence=SEQ, topology="linear", enzymes=None):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": topology,
            "enzymes": enzymes if enzymes is not None else [ECORI],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _partial(client, digest_id, cuts):
    return client.post(
        f"/digests/{digest_id}/partial-digest-distribution", json={"cuts": cuts}
    )


def test_single_cut_distribution(client):
    body = _create(client)
    resp = _partial(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 0.75}])
    assert resp.status_code == 200
    data = resp.json()

    assert data["id"] == body["id"]
    assert data["topology"] == "linear"
    assert data["cut_count"] == 1
    assert data["outcome_count"] == 2
    # 按概率降序：切开得 [5, 9]，未切开得全长 [14]
    assert data["outcomes"] == [
        {"fragment_lengths": [5, 9], "fragment_count": 2, "probability": 0.75},
        {"fragment_lengths": [14], "fragment_count": 1, "probability": 0.25},
    ]
    assert abs(data["total_probability"] - 1.0) <= 1e-9


def test_same_length_multiset_merges(client):
    # 只切 5 位得 [5, 9]，只切 9 位得 [9, 5]：长度多重集相同，合并累加
    body = _create(client, sequence=SEQ2, enzymes=[AVA, BVO])
    resp = _partial(
        client,
        body["id"],
        [
            {"top_cut": 5, "bottom_cut": 5, "probability": 0.8},
            {"top_cut": 9, "bottom_cut": 9, "probability": 0.5},
        ],
    )
    assert resp.status_code == 200
    data = resp.json()
    # 合并与排序：{5} 与 {9} 多重集相同合并为 [5,9]；两切口均成功得 [4,5,5]；
    # 均失败得全长 [14]。概率按降序排列
    expected = [
        ([5, 9], 2, 0.8 * 0.5 + 0.2 * 0.5),   # {5} 与 {9} 合并
        ([4, 5, 5], 3, 0.8 * 0.5),            # 两切口均成功
        ([14], 1, 0.2 * 0.5),                 # 均失败
    ]
    assert len(data["outcomes"]) == len(expected)
    for outcome, (lengths, count, prob) in zip(data["outcomes"], expected):
        assert outcome["fragment_lengths"] == lengths
        assert outcome["fragment_count"] == count
        assert outcome["probability"] == pytest.approx(prob, abs=1e-12)
    assert abs(data["total_probability"] - 1.0) <= 1e-9


def test_circular_merges_uncut_and_single_cut(client):
    # 环形：0 个或 1 个切口均为全长片段，三种组合合并；两切口得 [4, 12]
    body = _create(client, sequence=CIRC, topology="circular", enzymes=[ECORI, ECORV])
    resp = _partial(
        client,
        body["id"],
        [
            {"top_cut": 1, "bottom_cut": 5, "probability": 0.5},
            {"top_cut": 13, "bottom_cut": 13, "probability": 0.5},
        ],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["topology"] == "circular"
    assert data["outcomes"] == [
        {"fragment_lengths": [16], "fragment_count": 1, "probability": 0.75},
        {"fragment_lengths": [4, 12], "fragment_count": 2, "probability": 0.25},
    ]
    assert abs(data["total_probability"] - 1.0) <= 1e-9


def test_twelve_cuts_enumerated(client):
    # 13 个 A 上 pattern "AA"（offset 1/1）产生 12 个切口，恰为枚举上限
    poly_a = {"name": "PolyA", "pattern": "AA", "cut_top": 1, "cut_bottom": 1}
    body = _create(client, sequence="AAAAAAAAAAAAA", enzymes=[poly_a])
    assert len(body["cuts"]) == 12

    resp = _partial(
        client,
        body["id"],
        [{"top_cut": i, "bottom_cut": i, "probability": 0.5} for i in range(1, 13)],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["cut_count"] == 12
    # 4096 种组合按长度多重集合并后 = 13 的整数划分数
    assert data["outcome_count"] == 101
    assert abs(data["total_probability"] - 1.0) <= 1e-9
    # 降序排列
    probs = [o["probability"] for o in data["outcomes"]]
    assert probs == sorted(probs, reverse=True)


def test_too_many_cuts_returns_422(client):
    # 14 个 A -> 13 个切口，超过枚举上限
    poly_a = {"name": "PolyA", "pattern": "AA", "cut_top": 1, "cut_bottom": 1}
    body = _create(client, sequence="AAAAAAAAAAAAAA", enzymes=[poly_a])
    assert len(body["cuts"]) == 13
    resp = _partial(client, body["id"], [])
    assert resp.status_code == 422


def test_extreme_probabilities_filter_impossible_outcomes(client):
    body = _create(client)
    resp = _partial(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 1.0}])
    assert resp.status_code == 200
    data = resp.json()
    # 概率为 0 的组合不列入分布
    assert data["outcomes"] == [
        {"fragment_lengths": [5, 9], "fragment_count": 2, "probability": 1.0}
    ]
    assert data["total_probability"] == 1.0


def test_record_without_cuts(client):
    body = _create(client, enzymes=[ECORV])  # SEQ 上无 EcoRV 位点
    assert body["cuts"] == []
    resp = _partial(client, body["id"], [])
    assert resp.status_code == 200
    data = resp.json()
    assert data["cut_count"] == 0
    assert data["outcomes"] == [
        {"fragment_lengths": [14], "fragment_count": 1, "probability": 1.0}
    ]


def test_request_cut_order_does_not_matter(client):
    body = _create(client, sequence=SEQ2, enzymes=[AVA, BVO])
    cuts = [
        {"top_cut": 5, "bottom_cut": 5, "probability": 0.8},
        {"top_cut": 9, "bottom_cut": 9, "probability": 0.5},
    ]
    r1 = _partial(client, body["id"], cuts).json()
    r2 = _partial(client, body["id"], list(reversed(cuts))).json()
    assert r1 == r2


def test_partial_digest_is_read_only(client):
    body = _create(client, sequence=SEQ2, enzymes=[AVA, BVO])
    resp = _partial(
        client,
        body["id"],
        [
            {"top_cut": 5, "bottom_cut": 5, "probability": 0.8},
            {"top_cut": 9, "bottom_cut": 9, "probability": 0.5},
        ],
    )
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{body['id']}").json() == body
    assert len(client.get("/digests").json()) == 1


def test_missing_record_returns_404(client):
    resp = _partial(client, 999, [])
    assert resp.status_code == 404


def test_missing_coordinate_returns_422(client):
    body = _create(client, sequence=SEQ2, enzymes=[AVA, BVO])
    resp = _partial(client, body["id"], [{"top_cut": 5, "bottom_cut": 5, "probability": 0.5}])
    assert resp.status_code == 422
    assert "missing" in resp.json()["detail"]


def test_unknown_coordinate_returns_422(client):
    body = _create(client)
    resp = _partial(
        client,
        body["id"],
        [
            {"top_cut": 5, "bottom_cut": 9, "probability": 0.5},
            {"top_cut": 1, "bottom_cut": 1, "probability": 0.5},
        ],
    )
    assert resp.status_code == 422
    assert "unknown" in resp.json()["detail"]


def test_duplicate_coordinate_returns_422(client):
    body = _create(client)
    resp = _partial(
        client,
        body["id"],
        [
            {"top_cut": 5, "bottom_cut": 9, "probability": 0.5},
            {"top_cut": 5, "bottom_cut": 9, "probability": 0.6},
        ],
    )
    assert resp.status_code == 422
    assert "duplicate" in resp.json()["detail"]


def test_probability_out_of_range_returns_422(client):
    body = _create(client)
    for bad in (-0.1, 1.5):
        resp = _partial(
            client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": bad}]
        )
        assert resp.status_code == 422, (bad, resp.text)
