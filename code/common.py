"""Shared data access, plotting style and figure export."""

from pathlib import Path
import argparse

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "Figure_data.xlsx"
RESULTS = ROOT / "data" / "Results.xlsx"
OUTPUT = ROOT / "figures"

def frame(sheet):
    if not DATA.is_file():
        raise FileNotFoundError(f"Missing figure source workbook: {DATA}")
    return pd.read_excel(DATA, sheet_name=sheet, header=1, engine="openpyxl")


def table(component):
    """Read a numbered table component using the index in Results.xlsx."""
    with pd.ExcelFile(RESULTS, engine="openpyxl") as workbook:
        index = workbook.parse("Index").set_index("Table", verify_integrity=True)
        if component not in index.index:
            raise ValueError(f"No supplementary table component {component!r} in Results.xlsx")
        entry = index.loc[component]
        return workbook.parse(
            entry["Worksheet"],
            header=int(entry["Header row"]) - 1,
            nrows=int(entry["Data rows"]),
            usecols=range(int(entry["Columns"])),
        )

def save(fig, name, dpi=320):

    OUTPUT.mkdir(parents=True, exist_ok=True)

    for ext in ("png", "pdf", "svg"):

        fig.savefig(OUTPUT / f"{name}.{ext}", dpi=dpi)

    plt.close(fig)

def revision_style():

    mpl.rcdefaults()

    mpl.rcParams.update({

        "font.family": "DejaVu Sans", "font.size": 8,

        "axes.spines.top": False, "axes.spines.right": False,

        "pdf.fonttype": 42, "svg.fonttype": "none",

    })

def letter(ax, label):

    ax.text(-.14, 1.08, label, transform=ax.transAxes,

            fontweight="bold", fontsize=11)

def run(plot):
    global OUTPUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    OUTPUT = args.output
    plot()
