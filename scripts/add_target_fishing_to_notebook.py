#!/usr/bin/env python3
"""Embed the author-supplied target-fishing snapshot in the project Colab."""

import base64
import gzip
import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "docs/papers/terpene-docking-reproduction/notebooks/reproducibility.ipynb"
DOWNLOAD = ROOT / "docs/papers/terpene-docking-reproduction/downloads/reference-target-fishing.json"
SOURCE = Path("/Users/danielmcshan/GitHub/chat-terpedia-backend/research/reference-paper/derived/target-fishing-records.json")
TAG = "target-fishing-supplement"


def markdown(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {"tags": [TAG]}, "source": source.splitlines(keepends=True)}


def code(source: str) -> dict:
    return {
        "cell_type": "code", "execution_count": None, "metadata": {"tags": [TAG]},
        "outputs": [], "source": source.splitlines(keepends=True),
    }


def main() -> None:
    raw = SOURCE.read_bytes()
    expected = hashlib.sha256(raw).hexdigest()
    encoded = base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")
    DOWNLOAD.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE, DOWNLOAD)

    cells = [
        markdown("""## Author-supplied target-fishing reproduction

This section uses the checksummed Supplementary Information 1 and 2 snapshot supplied by the paper authors. It preserves 23 terpene identities, 1,546 row-level predictions from four tools, 260 author consensus rows, and the target-frequency table. The normalized snapshot is embedded so the notebook remains executable if an API or source spreadsheet changes.

The authors' consensus table remains authoritative. The independent reconstruction below groups targets by ChEMBL ID, then UniProt ID, then normalized target name and requires at least two distinct tools. Differences are reported rather than overwritten.
"""),
        code(f'''import gzip, unicodedata, re

EXPECTED_TARGET_FISHING_SHA256 = "{expected}"
TARGET_FISHING_BYTES = gzip.decompress(base64.b64decode("{encoded}"))
assert hashlib.sha256(TARGET_FISHING_BYTES).hexdigest() == EXPECTED_TARGET_FISHING_SHA256
TARGET_FISHING = json.loads(TARGET_FISHING_BYTES)

target_fishing_predictions = pd.DataFrame(TARGET_FISHING["predictions"])
author_consensus = pd.DataFrame(TARGET_FISHING["consensus_targets"])
target_frequencies = pd.DataFrame(TARGET_FISHING["target_frequencies"])

assert len(TARGET_FISHING["compounds"]) == 23
assert len(target_fishing_predictions) == 1546
assert set(target_fishing_predictions["tool"]) == {{"SwissTargetPrediction", "SEA", "STarFish", "PPB2"}}
assert target_frequencies.iloc[0]["target_name"] == "Cannabinoid receptor 2"
assert int(target_frequencies.iloc[0]["terpenes_predicting"]) == 18

print("Target-fishing snapshot:", EXPECTED_TARGET_FISHING_SHA256)
print(TARGET_FISHING["counts"])
'''),
        code('''def normalized_name(value):
    value = str(value or "").replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii").casefold()
    return re.sub(r"[^a-z0-9]", "", value)

def stable_target_key(row):
    return row.get("chembl_id") or row.get("uniprot_id") or normalized_name(row.get("target_name"))

target_fishing_predictions["ligand_key"] = target_fishing_predictions["terpene"].map(normalized_name)
target_fishing_predictions["target_key"] = target_fishing_predictions.apply(stable_target_key, axis=1)
author_consensus["ligand_key"] = author_consensus["terpene"].map(normalized_name)
author_consensus["target_key"] = author_consensus.apply(stable_target_key, axis=1)

def consensus_pairs(frame, omitted_tool=None, minimum_tools=2):
    current = frame if omitted_tool is None else frame[frame["tool"] != omitted_tool]
    grouped = (
        current.groupby(["ligand_key", "target_key"], as_index=False)["tool"]
        .agg(lambda values: sorted(set(values)))
        .rename(columns={"tool": "tools"})
    )
    grouped["tool_count"] = grouped["tools"].map(len)
    return grouped[grouped["tool_count"] >= minimum_tools].reset_index(drop=True)

recomputed_consensus = consensus_pairs(target_fishing_predictions)
author_pair_index = set(map(tuple, author_consensus[["ligand_key", "target_key"]].drop_duplicates().to_numpy()))
recomputed_pair_index = set(map(tuple, recomputed_consensus[["ligand_key", "target_key"]].to_numpy()))

consensus_comparison = {
    "recomputed_pairs": len(recomputed_pair_index),
    "author_pairs": len(author_pair_index),
    "overlap": len(recomputed_pair_index & author_pair_index),
    "precision_vs_author": len(recomputed_pair_index & author_pair_index) / len(recomputed_pair_index),
    "recall_vs_author": len(recomputed_pair_index & author_pair_index) / len(author_pair_index),
    "recomputed_only": len(recomputed_pair_index - author_pair_index),
    "author_only": len(author_pair_index - recomputed_pair_index),
}
pd.Series(consensus_comparison, name="value").to_frame()
'''),
        code('''TOOLS = ["SwissTargetPrediction", "SEA", "STarFish", "PPB2"]
ablation_rows = []
for omitted in TOOLS:
    reduced = consensus_pairs(target_fishing_predictions, omitted_tool=omitted)
    reduced_pairs = set(map(tuple, reduced[["ligand_key", "target_key"]].to_numpy()))
    ablation_rows.append({
        "omitted_tool": omitted,
        "consensus_pairs": len(reduced_pairs),
        "pairs_lost_from_recomputed_full": len(recomputed_pair_index - reduced_pairs),
        "author_pair_overlap": len(reduced_pairs & author_pair_index),
        "author_pair_recall": len(reduced_pairs & author_pair_index) / len(author_pair_index),
    })
leave_one_tool_out = pd.DataFrame(ablation_rows).sort_values("consensus_pairs", ascending=False)
leave_one_tool_out
'''),
        code('''top_five_targets = target_frequencies[target_frequencies["selected_top_five"]].copy()
display(top_five_targets[["rank", "target_name", "terpenes_predicting"]])

ax = top_five_targets.sort_values("terpenes_predicting").plot.barh(
    x="target_name", y="terpenes_predicting", legend=False, color="#28735d",
    title="Author-selected targets by predicted-terpene frequency",
)
ax.set_xlabel("Terpenes predicting target")
ax.set_ylabel("")
plt.tight_layout()
'''),
        code('''from pathlib import Path

target_fishing_output = Path("/content/terpedia-target-fishing" if "google.colab" in sys.modules else "target-fishing-outputs")
target_fishing_output.mkdir(parents=True, exist_ok=True)
target_fishing_predictions.to_csv(target_fishing_output / "normalized-predictions.csv", index=False)
author_consensus.to_csv(target_fishing_output / "author-consensus.csv", index=False)
recomputed_consensus.to_csv(target_fishing_output / "recomputed-consensus.csv", index=False)
leave_one_tool_out.to_csv(target_fishing_output / "leave-one-tool-out.csv", index=False)
top_five_targets.to_csv(target_fishing_output / "top-five-targets.csv", index=False)
(target_fishing_output / "consensus-comparison.json").write_text(json.dumps(consensus_comparison, indent=2))

print("Wrote", len(list(target_fishing_output.iterdir())), "target-fishing artifacts to", target_fishing_output)
'''),
    ]

    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    notebook["cells"] = [cell for cell in notebook["cells"] if TAG not in cell.get("metadata", {}).get("tags", [])]
    export_index = next(
        index for index, cell in enumerate(notebook["cells"])
        if cell["cell_type"] == "markdown" and "".join(cell.get("source", [])).strip().startswith("## Export")
    )
    notebook["cells"][export_index:export_index] = cells
    NOTEBOOK.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"notebook": str(NOTEBOOK), "cells_added": len(cells), "snapshot_sha256": expected}, indent=2))


if __name__ == "__main__":
    main()
