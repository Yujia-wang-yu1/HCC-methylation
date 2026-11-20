"""Write Table S2 to standard output as CSV.

Probe counts come from the HM450 processing summary in Supplementary Methods
S2. Sample counts come from the cohort data used in Figure 1.
Usage: python code/table_s2.py
"""

import csv
from pathlib import Path
import sys

from openpyxl import load_workbook


DATA = Path(__file__).resolve().parents[1] / "data"
SAMPLE_ROWS = [
    ("TCGA paired", "tcga_paired"),
    ("Maximum CNV/purity complete", "tcga_cnv_purity"),
    ("Broad DNA-proxy eligible", "tcga_dna"),
    ("Joint DNA/ABSOLUTE score complete", "tcga_dna_absolute"),
]


def main():
    with (DATA / "HM450_probe_counts.tsv").open(encoding="utf-8", newline="") as source:
        rows = [(row["measure"], int(row["n"])) for row in csv.DictReader(source, delimiter="\t")]
    workbook = load_workbook(DATA / "Figure_data.xlsx", read_only=True, data_only=True)
    try:
        sheet = workbook["Cohort framework"]
        counts = {row[0]: row[1] for row in sheet.iter_rows(min_row=3, values_only=True)}
        rows.extend((label, counts[key]) for label, key in SAMPLE_ROWS)
    finally:
        workbook.close()
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["measure", "n"])
    writer.writerows(rows)


if __name__ == "__main__":
    main()
