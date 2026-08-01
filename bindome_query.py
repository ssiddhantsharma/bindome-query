#!/usr/bin/env python3
"""bindome-query — pull designed-binder metrics for a list of targets from the
Human Bindome API (https://bindome.epfl.ch) and emit tidy CSVs + a figure.

The Human Bindome (Wenckstern et al., bioRxiv 2026.07.30.741542) is a
proteome-scale atlas of ~306k in-silico BindCraft binder candidates over ~8.3k
human proteins. For each target it exposes one record per designed binder with
AlphaFold-style confidence metrics and links to the predicted complex.

NOTE: these are computational design candidates, not experimentally validated
binders. There is no measured affinity here — treat the metrics as design-time
confidence only.

Usage:
    python bindome_query.py --targets targets.json --outdir out
    python bindome_query.py --targets targets.json --outdir out --top-structures 3
"""
from __future__ import annotations

import argparse
import difflib
import json
import re
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib.lines import Line2D

API = "https://bindome.epfl.ch"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

# (api_field, short label, higher_is_better)
METRICS = [
    ("average_i_pTM", "i_pTM", True),
    ("average_pTM", "pTM", True),
    ("average_i_pAE", "i_pAE", False),
    ("average_pAE", "pAE", False),
    ("average_pLDDT", "pLDDT", True),
    ("average_i_pLDDT", "i_pLDDT", True),
    ("average_ipSAE", "ipSAE", True),
]

# Metrics shown in the figure (the three most decision-relevant), in panel order.
FIGURE_METRICS = ["i_pTM", "i_pAE", "pLDDT"]

# BindCraft design names look like:
#   P24821_DOMAIN10-no_hotspot_l90_s989522_mpnn2_model2
ID_RE = re.compile(
    r"_(?P<domain>DOMAIN\d+)-(?P<hotspot>.+?)_l(?P<length>\d+)"
    r"_s(?P<seed>\d+)_mpnn(?P<mpnn>\d+)_model(?P<afmodel>\d+)$"
)

# ---- dataviz palette (validated reference palette, light mode) --------------
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
BOX_FILL = "#f0efec"
SURFACE = "#fcfcfb"
# Validated categorical palette (fixed order, never cycled). A 9th domain folds
# into "Other" (gray) rather than getting an unvalidated hue.
CATEGORICAL = ["#2a78d6", "#008300", "#e87ba4", "#eda100",
               "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"]
OTHER = "#898781"
THRESH_LINE = "#b5651d"  # dashed viability guide

# Typical BindCraft-style filter guides (normalized 0–1), per metric.
THRESHOLDS = {"i_pTM": 0.80, "i_pAE": 0.20, "pLDDT": 0.90}


# ---------------------------------------------------------------------------
# Networking
# ---------------------------------------------------------------------------
def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "bindome-query/1.0 (+https://github.com)"})
    return s


def _get(session: requests.Session, url: str, retries: int = 3, timeout: int = 45):
    """GET with retries. Returns None for any 4xx (absent/invalid resource) so one
    bad target never crashes the batch; retries only on network errors and 5xx."""
    last = None
    for i in range(retries):
        try:
            r = session.get(url, timeout=timeout)
        except requests.RequestException as e:  # transient network error
            last = e
            time.sleep(1.5 * (i + 1))
            continue
        if 400 <= r.status_code < 500:  # not found / bad accession — skip, don't crash
            return None
        if r.status_code >= 500:  # server error — retry
            last = requests.HTTPError(f"{r.status_code} Server Error for {url}")
            time.sleep(1.5 * (i + 1))
            continue
        return r
    raise RuntimeError(f"GET failed after {retries} tries: {url}\n  {last}")


def _norm(s: str) -> str:
    """'Claudin-6' -> 'claudin 6' (lowercase, non-alnum -> single space)."""
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _compact(s: str) -> str:
    """'PD-L1' -> 'pdl1' (alphanumerics only)."""
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def _candidate_names(gene_names: str, protein_name: str) -> list[str]:
    """All names an entry can be called by: gene symbols, the primary protein
    name, and every parenthesised synonym."""
    names: list[str] = []
    if gene_names:
        names += gene_names.split()
    if protein_name:
        names.append(re.split(r"[(\[]", protein_name)[0])  # primary, before '(' or '['
        names += re.findall(r"\(([^)]*)\)", protein_name)  # synonyms in parens
    return [n for n in (x.strip() for x in names) if n]


def pick_match(name: str, rows: list[dict]) -> dict:
    """Pure resolver logic (no network). `rows` are UniProt search hits with keys
    accession, id, gene_names, protein_name.

    Returns {'status', 'accession', 'display', 'candidates'} where status is:
      confident  — exactly one entry whose name/synonym/gene *exactly* matches
                   the input once normalized (spacing/case/punctuation ignored);
      ambiguous  — two or more distinct entries match exactly;
      notfound   — nothing matches exactly (candidates offered as suggestions).

    We never return a non-exact "closest" guess: fuzzy matching silently mistargets
    (e.g. 'claudn 6' is nearest to Claudin-3), which is unacceptable for design work.
    """
    want_n, want_c = _norm(name), _compact(name)
    exact: dict[str, dict] = {}
    scored: list[tuple[float, dict]] = []
    for row in rows:
        cnames = _candidate_names(row.get("gene_names", ""),
                                  row.get("protein_name", ""))
        display = cnames[0] if cnames else row.get("id", row["accession"])
        for nm in cnames:
            if _norm(nm) == want_n or (want_c and _compact(nm) == want_c):
                exact.setdefault(row["accession"],
                                 {"accession": row["accession"],
                                  "id": row.get("id"), "display": display})
                break
        # similarity to the input across all of this entry's names (for suggestions)
        sim = max((difflib.SequenceMatcher(None, want_n, _norm(nm)).ratio()
                   for nm in cnames), default=0.0)
        scored.append((sim, {"accession": row["accession"], "id": row.get("id"),
                             "display": display}))
    # suggestion list: most similar first, deduped by accession
    cand_view, seen_ac = [], set()
    for _, c in sorted(scored, key=lambda x: x[0], reverse=True):
        if c["accession"] not in seen_ac:
            seen_ac.add(c["accession"])
            cand_view.append(c)
    cand_view = cand_view[:6]
    if len(exact) == 1:
        hit = next(iter(exact.values()))
        return {"status": "confident", **hit, "candidates": cand_view}
    if len(exact) > 1:
        return {"status": "ambiguous", "accession": None, "display": None,
                "candidates": list(exact.values())}
    return {"status": "notfound", "accession": None, "display": None,
            "candidates": cand_view}


def _uniprot_search(session: requests.Session, query: str, size: int) -> list[dict]:
    fields = "accession,id,gene_names,protein_name"
    url = (f"{UNIPROT}/search?query={requests.utils.quote(query)}"
           f"&format=tsv&fields={fields}&size={size}")
    r = _get(session, url)
    if r is None:
        return []
    lines = [l for l in r.text.splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    keymap = {"Entry": "accession", "Entry Name": "id",
              "Gene Names": "gene_names", "Protein names": "protein_name"}
    out = []
    for line in lines[1:]:
        cells = line.split("\t")
        row = {keymap.get(h, h): (cells[i] if i < len(cells) else "")
               for i, h in enumerate(header)}
        if row.get("accession"):
            out.append(row)
    return out


def resolve_uniprot(session: requests.Session, name: str) -> dict:
    """Name -> reviewed human UniProt accession, tolerant of case / spacing /
    hyphen / gene-vs-name differences but never a fuzzy guess. Returns the dict
    from pick_match. A convenience only — prefer explicit `uniprot` in the
    targets file for anything load-bearing."""
    base = "AND organism_id:9606 AND reviewed:true"
    rows = _uniprot_search(session, f"{name} {base}", size=25)
    res = pick_match(name, rows)
    if res["status"] == "confident":
        return res
    # No exact hit from the plain query (often a typo). Widen with a wildcard on
    # the first alphabetic stem purely to offer "did you mean" suggestions.
    stem = re.findall(r"[A-Za-z]{3,}", name)
    if stem:
        wide = _uniprot_search(session, f"{stem[0][:4]}* {base}", size=100)
        # an exact match can still surface here (e.g. spacing-only mismatches)
        res2 = pick_match(name, wide)
        if res2["status"] == "confident":
            return res2
        if res["status"] == "notfound" and wide:
            res["candidates"] = res2["candidates"]
    return res


def fetch_target(session: requests.Session, ac: str) -> dict:
    """Return {'uniprot', 'id', 'records': [...]} for one accession."""
    r = _get(session, f"{API}/uniprot/{ac}/metrics")
    if r is None:
        return {"uniprot": ac, "id": None, "records": []}
    data = r.json() or {}
    entry = data.get("uniprot_entry") or {}
    return {
        "uniprot": entry.get("ac", ac),
        "id": entry.get("id"),
        "records": data.get("metrics") or [],
    }


# ---------------------------------------------------------------------------
# Parsing / tables
# ---------------------------------------------------------------------------
def parse_design_name(model_identifier: str) -> dict:
    m = ID_RE.search(model_identifier or "")
    if not m:
        return {k: None for k in ("domain", "hotspot", "binder_length", "seed",
                                  "mpnn_variant", "af_model")}
    g = m.groupdict()
    return {
        "domain": g["domain"],
        "hotspot": g["hotspot"],
        "binder_length": int(g["length"]),
        "seed": int(g["seed"]),
        "mpnn_variant": int(g["mpnn"]),
        "af_model": int(g["afmodel"]),
    }


def designs_frame(targets: list[dict]) -> pd.DataFrame:
    rows = []
    for t in targets:
        for rec in t["records"]:
            row = {
                "target_name": t["name"],
                "uniprot": t["uniprot"],
                "uniprot_id": t["id"],
                "model_identifier": rec.get("model_identifier"),
                "uniprot_start": rec.get("uniprot_start"),
                "uniprot_end": rec.get("uniprot_end"),
                "coverage": rec.get("coverage"),
            }
            row.update(parse_design_name(rec.get("model_identifier", "")))
            for api_field, short, _ in METRICS:
                row[short] = rec.get(api_field)
            row["model_url"] = rec.get("model_url")
            row["pae_url"] = rec.get("pae_url")
            rows.append(row)
    cols = [
        "target_name", "uniprot", "uniprot_id", "model_identifier",
        "domain", "hotspot", "binder_length", "seed", "mpnn_variant", "af_model",
        "uniprot_start", "uniprot_end", "coverage",
        *[s for _, s, _ in METRICS], "model_url", "pae_url",
    ]
    return pd.DataFrame(rows, columns=cols)


def summary_frame(targets: list[dict], designs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for t in targets:
        d = designs[designs["uniprot"] == t["uniprot"]]
        row = {
            "target_name": t["name"],
            "uniprot": t["uniprot"],
            "uniprot_id": t["id"],
            "n_binders": len(d),
            "n_regions": d[["uniprot_start", "uniprot_end"]].drop_duplicates().shape[0]
            if len(d) else 0,
        }
        for _, short, higher in METRICS:
            vals = d[short].dropna() if len(d) else pd.Series(dtype=float)
            if len(vals):
                row[f"best_{short}"] = vals.max() if higher else vals.min()
                row[f"median_{short}"] = vals.median()
            else:
                row[f"best_{short}"] = np.nan
                row[f"median_{short}"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("best_i_pTM", ascending=False,
                                          na_position="last").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def _domain_colors(have: pd.DataFrame) -> dict:
    """Map each target domain to a fixed categorical hue (most frequent first);
    domains past the palette fold into a single gray 'Other'."""
    doms = (have["domain"].fillna("?").value_counts().index.tolist())
    mapping = {}
    for i, d in enumerate(doms):
        mapping[d] = CATEGORICAL[i] if i < len(CATEGORICAL) else OTHER
    return mapping


def make_figure(designs: pd.DataFrame, empty_targets: list[str], out_png: Path,
                updated: str) -> None:
    have = designs.dropna(subset=["i_pTM"]).copy()
    if have.empty:
        print("  [figure] no binders across any target — skipping figure")
        return
    have["domain"] = have["domain"].fillna("?")

    order = (have.groupby("target_name")["i_pTM"].median()
             .sort_values().index.tolist())
    counts = have.groupby("target_name")["i_pTM"].size()
    ylabels = [f"{name}  (n={counts[name]})" for name in order]
    ypos = {name: i for i, name in enumerate(order)}

    dom_color = _domain_colors(have)
    multi_domain = len(dom_color) > 1

    metrics = [(s, h) for _, s, h in METRICS if s in FIGURE_METRICS]
    arrows = {True: "↑ higher better", False: "↓ lower better"}

    fig, axes = plt.subplots(
        1, len(metrics), sharey=True,
        figsize=(3.6 * len(metrics) + 1.2, 0.52 * len(order) + 2.6),
    )
    if len(metrics) == 1:
        axes = [axes]
    fig.patch.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)

    for ax, (short, higher) in zip(axes, metrics):
        ax.set_facecolor(SURFACE)
        thr = THRESHOLDS.get(short)
        if thr is not None:
            ax.axvline(thr, color=THRESH_LINE, ls=(0, (4, 3)), lw=1.1,
                       zorder=2, alpha=0.7)
        for name in order:
            sub = have.loc[have["target_name"] == name]
            vals = sub[short].dropna().values
            y = ypos[name]
            if not len(vals):
                continue
            ax.boxplot(
                vals, positions=[y], orientation="horizontal", widths=0.55,
                patch_artist=True, showfliers=False, zorder=1,
                medianprops=dict(color=INK2, lw=1.4),
                whiskerprops=dict(color=BASE, lw=1.0),
                capprops=dict(color=BASE, lw=1.0),
                boxprops=dict(facecolor=BOX_FILL, edgecolor=BASE, lw=1.0),
            )
            jit = rng.uniform(-0.16, 0.16, size=len(vals))
            colors = [dom_color[d] for d in sub["domain"].values]
            ax.scatter(vals, np.full(len(vals), y) + jit, s=24, c=colors,
                       edgecolor="white", linewidth=0.3, zorder=3, alpha=0.95)
            # median value label above the box
            med = float(np.median(vals))
            ax.text(med, y + 0.34, f"{med:.2f}", fontsize=7.5, color=INK2,
                    ha="center", va="bottom", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.1", fc=SURFACE, ec="none",
                              alpha=0.7))
        ax.set_title(f"{short}\n{arrows[higher]}", fontsize=11, color=INK,
                     pad=8, linespacing=1.4)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(BASE)

    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(ylabels, fontsize=9, color=INK)
    axes[0].set_ylim(-0.7, len(order) + 0.1)

    fig.suptitle("Human Bindome — designed-binder confidence by target",
                 fontsize=14, color=INK, x=0.02, ha="left", y=0.99)

    # domain legend (only when there's more than one domain to distinguish),
    # ordered naturally by domain number rather than by frequency
    if multi_domain:
        def _dnum(d):
            m = re.search(r"\d+", d)
            return int(m.group()) if m else 999
        items = sorted(dom_color.items(), key=lambda kv: _dnum(kv[0]))
        handles = [Line2D([0], [0], marker="o", ls="", markersize=7,
                          markerfacecolor=c, markeredgecolor="white",
                          markeredgewidth=0.4, label=d)
                   for d, c in items]
        fig.legend(handles=handles, title="target domain", loc="upper left",
                   bbox_to_anchor=(0.01, 0.955), ncol=min(len(handles), 8),
                   frameon=False, fontsize=8, title_fontsize=8,
                   handletextpad=0.3, columnspacing=1.1)

    guide = ", ".join(f"{k}{'≥' if h else '≤'}{THRESHOLDS[k]:g}"
                      for _, k, h in METRICS if k in THRESHOLDS)
    lines = [
        f"Source: bindome.epfl.ch  ·  metrics normalized 0–1 as returned by the "
        f"API  ·  in-silico BindCraft candidates (no measured affinity)  ·  data {updated}",
        f"dashed line = typical filter: {guide}",
    ]
    if empty_targets:
        lines.append(f"no Bindome binders for: {', '.join(empty_targets)}")
    fig.text(0.02, 0.005, "\n".join(lines), fontsize=8, color=MUTED,
             ha="left", va="bottom", linespacing=1.5)

    top = 0.90 if multi_domain else 0.93
    fig.tight_layout(rect=[0, 0.06 + 0.02 * len(lines), 1, top])
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"  [figure] wrote {out_png}")


# ---------------------------------------------------------------------------
# Structures (optional)
# ---------------------------------------------------------------------------
def download_top_structures(session: requests.Session, designs: pd.DataFrame,
                            outdir: Path, top_n: int) -> None:
    sdir = outdir / "structures"
    sdir.mkdir(parents=True, exist_ok=True)
    got = 0
    for name, d in designs.dropna(subset=["i_pTM"]).groupby("target_name"):
        top = d.sort_values("i_pTM", ascending=False).head(top_n)
        for _, row in top.iterrows():
            url = row["model_url"]
            if not isinstance(url, str):
                continue
            dest = sdir / f"{row['uniprot']}__{Path(url).name}"
            r = _get(session, url)
            if r is None:
                print(f"  [struct] MISSING {url}")
                continue
            dest.write_bytes(r.content)
            got += 1
    print(f"  [struct] downloaded {got} top structure(s) into {sdir}")


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------
def load_targets(path: Path) -> list[dict]:
    raw = json.loads(path.read_text())
    items = raw["targets"] if isinstance(raw, dict) and "targets" in raw else raw
    out = []
    for it in items:
        if isinstance(it, str):
            out.append({"name": it, "uniprot": None})
        else:
            out.append({"name": it.get("name") or it.get("uniprot"),
                        "uniprot": it.get("uniprot")})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", required=True, type=Path,
                    help="JSON file: list of {name, uniprot} (see targets.example.json)")
    ap.add_argument("--outdir", default=Path("out"), type=Path)
    ap.add_argument("--top-structures", type=int, default=0,
                    help="also download the top-N binder CIFs per target (0 = off)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()
    targets = load_targets(args.targets)

    resolved: list[dict] = []
    seen: set[str] = set()
    for t in targets:
        ac = t["uniprot"]
        if not ac:
            res = resolve_uniprot(session, t["name"])
            if res["status"] == "confident":
                ac = res["accession"]
                print(f"  [resolve] '{t['name']}' -> {ac} {res['display']} (verify)")
            else:
                sugg = ", ".join(f"{c['display']} ({c['accession']})"
                                 for c in res["candidates"][:5]) or "none found"
                why = ("ambiguous — matches several entries"
                       if res["status"] == "ambiguous" else "no exact match")
                print(f"  [resolve] '{t['name']}': {why}; did you mean: {sugg}")
                print(f"            -> skipping (add an explicit \"uniprot\" to use one)")
                continue
        if ac in seen:  # same accession twice would double counts downstream
            print(f"  [skip]  {t['name']:16} {ac:10} duplicate accession — skipping")
            continue
        seen.add(ac)
        got = fetch_target(session, ac)
        got["name"] = t["name"]
        n = len(got["records"])
        flag = "" if n else "  (no binders in Bindome)"
        print(f"  [fetch] {t['name']:16} {ac:10} {n:>4} binders{flag}")
        resolved.append(got)

    if not resolved:
        print("No targets resolved. Nothing to do.")
        return 1

    designs = designs_frame(resolved)
    summary = summary_frame(resolved, designs)
    updated = time.strftime("%Y-%m-%d")

    designs.to_csv(args.outdir / "bindome_designs.csv", index=False)
    summary.to_csv(args.outdir / "bindome_summary.csv", index=False)
    print(f"\n  [csv] {len(designs)} designs -> {args.outdir/'bindome_designs.csv'}")
    print(f"  [csv] {len(summary)} targets -> {args.outdir/'bindome_summary.csv'}")

    empty = [t["name"] for t in resolved if not t["records"]]
    make_figure(designs, empty, args.outdir / "bindome_metrics.png", updated)

    if args.top_structures > 0:
        download_top_structures(session, designs, args.outdir, args.top_structures)

    # Console summary
    print("\nSummary (best i_pTM shown):")
    view = summary[["target_name", "uniprot", "n_binders", "n_regions",
                    "best_i_pTM", "median_i_pTM", "best_i_pAE"]]
    print(view.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
