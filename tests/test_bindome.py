"""Unit tests for bindome-query. No network: every test exercises pure logic or
mocked inputs, so `pytest` runs green in CI without hitting the live API."""
import pandas as pd
import pytest

import bindome_query as bq


# --------------------------------------------------------------------------
# design-name parsing
# --------------------------------------------------------------------------
def test_parse_design_name_valid():
    got = bq.parse_design_name(
        "P24821_DOMAIN10-no_hotspot_l90_s989522_mpnn2_model2")
    assert got == {
        "domain": "DOMAIN10", "hotspot": "no_hotspot", "binder_length": 90,
        "seed": 989522, "mpnn_variant": 2, "af_model": 2,
    }


def test_parse_design_name_malformed_is_nulls_not_crash():
    got = bq.parse_design_name("totally-not-a-bindcraft-name")
    assert set(got) == {"domain", "hotspot", "binder_length", "seed",
                        "mpnn_variant", "af_model"}
    assert all(v is None for v in got.values())


def test_parse_design_name_empty():
    assert bq.parse_design_name("")["domain"] is None


# --------------------------------------------------------------------------
# normalization helpers
# --------------------------------------------------------------------------
@pytest.mark.parametrize("raw,norm,compact", [
    ("Claudin-6", "claudin 6", "claudin6"),
    ("PD-L1", "pd l1", "pdl1"),
    ("  HER2 ", "her2", "her2"),
])
def test_norm_and_compact(raw, norm, compact):
    assert bq._norm(raw) == norm
    assert bq._compact(raw) == compact


def test_candidate_names_extracts_gene_primary_and_synonyms():
    names = bq._candidate_names(
        "CD274 PDCD1LG1 PDL1",
        "Programmed cell death 1 ligand 1 (PD-L1) (B7-H1) [Cleaved]")
    assert "CD274" in names and "PDL1" in names
    assert "Programmed cell death 1 ligand 1" in names
    assert "PD-L1" in names and "B7-H1" in names


# --------------------------------------------------------------------------
# resolver (pure) — the safety-critical part
# --------------------------------------------------------------------------
CLDN6 = {"accession": "P56747", "id": "CLD6_HUMAN",
         "gene_names": "CLDN6", "protein_name": "Claudin-6 (Skullin)"}
CLDN9 = {"accession": "O95484", "id": "CLD9_HUMAN",
         "gene_names": "CLDN9", "protein_name": "Claudin-9"}
CLDN3 = {"accession": "O15551", "id": "CLD3_HUMAN",
         "gene_names": "CLDN3", "protein_name": "Claudin-3"}
PDL1 = {"accession": "Q9NZQ7", "id": "PD1L1_HUMAN",
        "gene_names": "CD274 PDCD1LG1 PDL1",
        "protein_name": "Programmed cell death 1 ligand 1 (PD-L1) (B7-H1)"}


def test_pick_match_prefers_exact_over_relevance():
    # 'claudin 6' must resolve to Claudin-6 even if Claudin-9 is listed first.
    res = bq.pick_match("claudin 6", [CLDN9, CLDN6])
    assert res["status"] == "confident"
    assert res["accession"] == "P56747"


def test_pick_match_synonym_and_compact():
    assert bq.pick_match("PD-L1", [PDL1])["accession"] == "Q9NZQ7"   # synonym
    assert bq.pick_match("pdl1", [PDL1])["accession"] == "Q9NZQ7"    # gene, compact


def test_pick_match_gene_symbol():
    assert bq.pick_match("CLDN6", [CLDN6, CLDN9])["accession"] == "P56747"


def test_pick_match_ambiguous():
    dup = {"accession": "X1", "id": "X1", "gene_names": "FOO",
           "protein_name": "Foo"}
    dup2 = {"accession": "X2", "id": "X2", "gene_names": "BAR",
            "protein_name": "Foo"}
    res = bq.pick_match("foo", [dup, dup2])
    assert res["status"] == "ambiguous"
    assert res["accession"] is None
    assert {c["accession"] for c in res["candidates"]} == {"X1", "X2"}


def test_pick_match_typo_never_guesses_but_suggests_closest_first():
    # 'claudn 6' has NO exact match -> notfound, but Claudin-6 is the top hint.
    res = bq.pick_match("claudn 6", [CLDN3, CLDN6, CLDN9])
    assert res["status"] == "notfound"
    assert res["accession"] is None
    assert res["candidates"][0]["accession"] == "P56747"


def test_pick_match_empty_rows():
    res = bq.pick_match("anything", [])
    assert res["status"] == "notfound"
    assert res["candidates"] == []


# --------------------------------------------------------------------------
# tables
# --------------------------------------------------------------------------
def _target(name, ac, n, ipae_lo=0.2, iptm_hi=0.8):
    """Build a fetched-target dict with n fake binder records."""
    recs = []
    for i in range(n):
        recs.append({
            "model_identifier": f"{ac}_DOMAIN1-no_hotspot_l90_s{i}_mpnn2_model1",
            "uniprot_start": 10, "uniprot_end": 100, "coverage": 0.3,
            "average_i_pTM": iptm_hi - 0.01 * i, "average_pTM": 0.75,
            "average_i_pAE": ipae_lo + 0.01 * i, "average_pAE": 0.2,
            "average_pLDDT": 0.9, "average_i_pLDDT": 0.5, "average_ipSAE": 0.9,
            "model_url": f"https://x/{ac}_{i}.cif", "pae_url": f"https://x/{ac}_{i}.json",
        })
    return {"name": name, "uniprot": ac, "id": f"{name}_ID", "records": recs}


def test_designs_frame_columns_and_parse():
    d = bq.designs_frame([_target("A", "P1", 2)])
    assert len(d) == 2
    for col in ("target_name", "uniprot", "domain", "binder_length",
                "i_pTM", "i_pAE", "model_url"):
        assert col in d.columns
    assert d["domain"].iloc[0] == "DOMAIN1"
    assert d["binder_length"].iloc[0] == 90


def test_summary_best_uses_min_for_lower_is_better():
    t = _target("A", "P1", 3, ipae_lo=0.2, iptm_hi=0.9)  # i_pAE 0.20,0.21,0.22
    d = bq.designs_frame([t])
    s = bq.summary_frame([t], d)
    row = s[s["uniprot"] == "P1"].iloc[0]
    assert row["best_i_pTM"] == pytest.approx(0.90)   # higher-is-better -> max
    assert row["best_i_pAE"] == pytest.approx(0.20)   # lower-is-better  -> min
    assert row["n_binders"] == 3


def test_summary_zero_binder_target_surfaced_not_dropped():
    t = {"name": "Z", "uniprot": "P0", "id": None, "records": []}
    d = bq.designs_frame([t])
    s = bq.summary_frame([t], d)
    assert (s["uniprot"] == "P0").any()
    row = s[s["uniprot"] == "P0"].iloc[0]
    assert row["n_binders"] == 0
    assert pd.isna(row["best_i_pTM"])


def test_summary_no_double_count_for_same_accession():
    # Two fetched entries with the SAME accession must not inflate each other.
    t = _target("A", "P1", 4)
    d = bq.designs_frame([t])          # a single fetch, 4 records
    s = bq.summary_frame([t], d)
    assert int(s[s["uniprot"] == "P1"]["n_binders"].iloc[0]) == 4


# --------------------------------------------------------------------------
# input loading
# --------------------------------------------------------------------------
def test_load_targets_variants(tmp_path):
    import json
    p = tmp_path / "t.json"
    p.write_text(json.dumps({"targets": [
        {"name": "A", "uniprot": "P1"},
        {"uniprot": "P2"},
        "P3",
    ]}))
    got = bq.load_targets(p)
    assert got[0] == {"name": "A", "uniprot": "P1"}
    assert got[1] == {"name": "P2", "uniprot": "P2"}   # name falls back to AC
    assert got[2] == {"name": "P3", "uniprot": None}   # bare string -> resolve later


def test_load_targets_bare_list(tmp_path):
    import json
    p = tmp_path / "t.json"
    p.write_text(json.dumps([{"name": "A", "uniprot": "P1"}]))
    assert bq.load_targets(p)[0]["uniprot"] == "P1"


# --------------------------------------------------------------------------
# figure smoke test (Agg backend, writes a file)
# --------------------------------------------------------------------------
def test_make_figure_writes_png(tmp_path):
    d = bq.designs_frame([_target("A", "P1", 5), _target("B", "P2", 6)])
    out = tmp_path / "fig.png"
    bq.make_figure(d, empty_targets=["EGFR"], out_png=out, updated="2026-01-01")
    assert out.exists() and out.stat().st_size > 1000


def test_make_figure_no_binders_skips_cleanly(tmp_path):
    d = bq.designs_frame([{"name": "Z", "uniprot": "P0", "id": None, "records": []}])
    out = tmp_path / "fig.png"
    bq.make_figure(d, empty_targets=["Z"], out_png=out, updated="2026-01-01")
    assert not out.exists()   # nothing to plot -> no file, no crash
