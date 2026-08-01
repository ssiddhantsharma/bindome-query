# bindome-query

[![CI](https://github.com/ssiddhantsharma/bindome-query/actions/workflows/ci.yml/badge.svg)](https://github.com/ssiddhantsharma/bindome-query/actions/workflows/ci.yml)

Query the [Human Bindome](https://bindome.epfl.ch) for a list of protein targets.
It pulls the metrics for every designed binder, writes a per-target summary, and
draws a figure of how those designs score.

The Human Bindome (Wenckstern *et al.*, bioRxiv
[2026.07.30.741542](https://www.biorxiv.org/content/10.64898/2026.07.30.741542v1))
is a proteome-scale atlas of about 306,000 in-silico
[BindCraft](https://github.com/martinpacesa/BindCraft) binder candidates covering
roughly 8,300 human proteins (40.9% of the proteome). Each target comes with one
record per designed binder, and each record carries AlphaFold-style confidence
metrics plus links to the predicted complex.

> Heads up: these are computational design candidates, not experimentally
> validated binders. There is no measured affinity in the Bindome. The metrics
> (`i_pTM`, `i_pAE`, `pLDDT`, and so on) are design-time confidence only, so use
> them to triage, not to conclude that a target "has a binder".

![example figure](docs/bindome_metrics.png)

## Install

```bash
pip install -r requirements.txt
```

## Usage

```bash
python bindome_query.py --targets targets.example.json --outdir out
```

Add `--top-structures N` to also download the top N binder CIFs per target
(ranked by `i_pTM`) into `out/structures/`.

### Targets file

A JSON list of `{name, uniprot}`. A bare list of accession strings works too.

```json
{
  "targets": [
    { "name": "PD-L1",  "uniprot": "Q9NZQ7" },
    { "name": "VEGF-A", "uniprot": "P15692" },
    { "name": "IL-2",   "uniprot": "P60568" }
  ]
}
```

### Name resolution

If you leave out `uniprot`, the name is looked up against UniProt (reviewed,
human). The lookup is tolerant of case, spacing, hyphens, and gene vs protein
names, so `CLDN6`, `Claudin-6`, `claudin 6`, `PD-L1`, `pdl1`, and `HER2` all
resolve correctly. It matches on the exact normalized name, so it never silently
picks the wrong protein. (Naive fuzzy search maps `claudin 6` to Claudin-9 and
`her-2` to HERC3; this tool does not.)

For a real typo it will not guess. It prints the closest "did you mean"
candidates and skips that target:

```
[resolve] 'claudn 6': no exact match; did you mean: CLDN6 (P56747), CLDN16 (Q9Y5I7), ...
          -> skipping (add an explicit "uniprot" to use one)
```

UniProt has no fuzzy search, and for design work picking the wrong target is
worse than a clean miss. If a target matters, put its accession in the file.
That is the key the Bindome API uses anyway.

## Outputs

| File | Contents |
|---|---|
| `out/bindome_designs.csv` | one row per designed binder: target, the UniProt region hit (`uniprot_start/end`), parsed design fields (`domain`, `binder_length`, `seed`, `mpnn_variant`, `af_model`), all seven metrics, and `model_url`/`pae_url` |
| `out/bindome_summary.csv` | one row per target: `n_binders`, `n_regions`, and the best and median of every metric |
| `out/bindome_metrics.png` | two metric-vs-metric gate scatters (`i_pTM` vs `i_pAE`, `pLDDT` vs `ipSAE`) colored by target, with the passing corner shaded and dashed cutoff lines, plus a target key and a per-target pass-rate table |

## Metrics

All values are normalized to 0 to 1 as returned by the API. Higher is better for
`i_pTM`, `pTM`, `pLDDT`, `i_pLDDT`, and `ipSAE`; lower is better for `i_pAE` and
`pAE`. The `uniprot_start` and `uniprot_end` columns tell you which region of the
target each binder was designed against, which is handy for checking whether a
candidate overlaps an epitope you care about.

Cutoffs used in the figure (normalized units): `i_pTM ≥ 0.85`, `i_pAE ≤ 0.24`,
`pLDDT ≥ 0.90`, `ipSAE ≥ 0.90`. One thing to watch: `i_pAE` is normalized to
about Å/25 (I calibrated this against the raw PAE matrices), so `i_pAE ≤ 0.24` is
roughly 6 Å, whereas `0.6` would be around 15 Å. To change any of them, edit
`THRESHOLDS` in `bindome_query.py`. And again, these are in-silico BindCraft
design metrics, not measured affinity.

## API notes

The tool uses one unauthenticated endpoint:

```
GET https://bindome.epfl.ch/uniprot/{UNIPROT_AC}/metrics
```

There is no batch or list endpoint, so the tool just loops over your targets. If
you want everything in bulk (all sequences, structures, and ML splits), the
authors publish the full dataset on Hugging Face:
[`wjulius/HumanBindome`](https://huggingface.co/datasets/wjulius/HumanBindome).

## Tests

```bash
pip install -r requirements-dev.txt
pytest -q
```

The tests run fully offline, with no API calls. They cover design-name parsing,
the name-resolver safety logic, the metric tables (including the no-double-count
and lower-is-better handling), and a figure smoke test. CI runs them on Python
3.10 to 3.12 (`.github/workflows/ci.yml`).

## License

MIT. The Human Bindome data is CC-BY-4.0, so please cite the paper if you use it.
