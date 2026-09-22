"""部分酶切（partial digest）概率分布分析。

对已保存记录的每个物理切口赋予独立的切割概率，枚举全部"切开/未切开"
组合（至多 12 个切口，即 4096 种），复用线性/环形片段重建逻辑得到每种
组合的片段长度多重集；多重集相同的结果合并、概率累加，按概率降序返回。
整个分析只在内存中计算，不读写数据库，原记录保持不变。
"""

from __future__ import annotations

from math import fsum

from .digest import build_fragments

# 枚举的切口数上限：2^12 = 4096 种组合，计算量可控
MAX_PARTIAL_CUTS = 12


def align_cut_probabilities(saved_cuts, requested):
    """校验请求坐标与已保存物理切口完全一致，并返回对齐后的概率列表。

    saved_cuts 为已保存切口（按 (top_cut, bottom_cut) 排序，物理坐标唯一）；
    requested 为请求中的 CutProbabilityIn 列表。校验通过时返回
    (按 saved_cuts 顺序排列的概率列表, None)；否则返回 (None, 错误描述)，
    由路由转为 422。
    """
    seen: dict[tuple[int, int], float] = {}
    duplicates: list[tuple[int, int]] = []
    for item in requested:
        key = (item.top_cut, item.bottom_cut)
        if key in seen:
            duplicates.append(key)
        seen[key] = item.probability
    if duplicates:
        return None, "duplicate cut coordinates: " + ", ".join(
            f"(top_cut={t}, bottom_cut={b})" for t, b in sorted(set(duplicates))
        )

    saved_keys = [(c["top_cut"], c["bottom_cut"]) for c in saved_cuts]
    saved_set = set(saved_keys)
    missing = [k for k in saved_keys if k not in seen]
    unknown = [k for k in seen if k not in saved_set]
    if missing or unknown:
        parts = []
        if missing:
            parts.append(
                "missing "
                + ", ".join(f"(top_cut={t}, bottom_cut={b})" for t, b in missing)
            )
        if unknown:
            parts.append(
                "unknown "
                + ", ".join(f"(top_cut={t}, bottom_cut={b})" for t, b in unknown)
            )
        return None, "cut coordinates do not match the saved record: " + "; ".join(parts)

    return [seen[k] for k in saved_keys], None


def partial_digest_distribution(sequence, topology, cuts, probabilities):
    """枚举全部切开/未切开组合，按片段长度多重集合并概率。

    cuts 为已保存切口（按 (top_cut, bottom_cut) 排序），probabilities
    与其一一对应。返回 {"outcomes": [...], "total_probability": float}：
    outcomes 按概率降序、概率相同按长度多重集升序排列（输出完全确定）；
    概率为 0 的组合（某切口概率取 0 或 1 时）不列入分布。
    """
    n = len(cuts)
    grouped: dict[tuple[int, ...], list[float]] = {}
    for mask in range(1 << n):
        prob = 1.0
        active = []
        for i, cut in enumerate(cuts):
            if mask >> i & 1:
                prob *= probabilities[i]
                active.append(cut)
            else:
                prob *= 1.0 - probabilities[i]
        if prob == 0.0:
            continue  # 不可能出现的组合
        fragments = build_fragments(sequence, topology, active)
        key = tuple(sorted(f["length"] for f in fragments))
        grouped.setdefault(key, []).append(prob)

    outcomes = [
        {
            "fragment_lengths": list(lengths),
            "fragment_count": len(lengths),
            "probability": fsum(probs),
        }
        for lengths, probs in grouped.items()
    ]
    # 概率降序为主键、长度多重集升序为次键，保证相同输入的输出必然一致
    outcomes.sort(key=lambda o: (-o["probability"], o["fragment_lengths"]))
    return {
        "outcomes": outcomes,
        "total_probability": fsum(o["probability"] for o in outcomes),
    }
