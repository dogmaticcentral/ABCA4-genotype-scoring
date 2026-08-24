#!/usr/bin/env python3
"""Port of the original abca4_60_production/generic/33_onset_age_vs_score.py to SQLite.

Correlate genotype functionality scores with the reported age of onset for the
cases of one publication, print the Spearman/Pearson correlations, and plot
onset age vs score (scatter, violin, or whisker).

The scores must already be stored in the database - run score_n_store.py first
(and with --assume_dosage_compensation if you want -s score_w_dosage_compensation).

Unlike the scoring pipeline itself, this script needs numpy, scipy and matplotlib.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, spearmanr, pearsonr

from abca4_score.abca4_queries import genotype_has_bad_allele
from abca4_score.sqlite_utils import connect, error_intolerant_search, hard_landing_search

DEFAULT_DB = Path(__file__).parent / "data" / "abca4_pub70_test.db"


# ---------------------------------------------------------------------------
@dataclass
class DensityEstimate:
    """Result of the KDE estimation."""

    grid: np.ndarray
    density: np.ndarray
    boundaries: np.ndarray
    peaks: np.ndarray


# ---------------------------------------------------------------------------
def estimate_density(x: np.ndarray, bandwidth: float | None = None, grid_size: int = 1000,
                     prominence: float = 0.03) -> DensityEstimate:
    """
    Estimate the KDE and locate cluster boundaries.

    Parameters
    ----------
    bandwidth
        KDE bandwidth multiplier (None = Scott).
    prominence
        Minimum peak prominence relative to the maximum density.
    """

    kde = gaussian_kde(x, bw_method=bandwidth)

    xmin = x.min()
    xmax = x.max()

    padding = 0.05 * (xmax - xmin)
    grid = np.linspace(xmin - padding, xmax + padding, grid_size)
    density = kde(grid)

    peaks, _ = find_peaks(density, prominence=prominence * density.max())
    boundaries = []

    for left, right in zip(peaks[:-1], peaks[1:]):
        minimum = left + np.argmin(density[left:right + 1])
        boundaries.append(grid[minimum])

    return DensityEstimate(grid=grid, density=density, peaks=peaks, boundaries=np.asarray(boundaries))


def _bin_onset_ages_by_score(scores, onset_ages):
    """Group onset ages into score bins of width 0.1 from 0 to 1."""
    bin_edges = np.arange(0, 1.01, 0.1)
    bin_labels = [f"{bin_edges[i]:.2f}-{bin_edges[i + 1]:.2f}"
                  for i in range(len(bin_edges) - 1)]

    # Digitize scores into bins (1-based indexing)
    bin_indices = np.digitize(scores, bin_edges, right=False)

    binned_data = []
    positions = []
    for i in range(1, len(bin_edges)):  # bins 1 to 20
        mask = (bin_indices == i)
        if np.any(mask):
            binned_data.append(onset_ages[mask])
            positions.append(i)

    return binned_data, positions, bin_labels


def _plot_boxplot_or_violin(ax: Axes, scores, onset_ages, plot_type: str):
    """Plot onset ages grouped into score bins as a boxplot or violin plot."""
    binned_data, positions, bin_labels = _bin_onset_ages_by_score(scores, onset_ages)

    fig, ax = plt.subplots(figsize=(8, 6))

    if plot_type == "violin-basic":
        ax.violinplot(binned_data, positions=positions, showmedians=True, showextrema=True)
    else:
        ax.boxplot(binned_data,
                   positions=positions,
                   patch_artist=True,
                   whiskerprops=dict(color='black', linewidth=1.5),
                   flierprops=dict(marker='o', markersize=4, alpha=0.6),
                   medianprops=dict(color='red', linewidth=2))

    ax.set_xticks(range(1, len(bin_labels) + 1))
    ax.set_xticklabels(bin_labels, rotation=45, ha='right')
    ax.set_xlim(0, len(positions) + 1)

    ax.set_xlabel('Score Bins')
    ax.set_ylabel('Onset Age')
    ax.set_title('Distribution of Onset Ages by Score Bins')

    ax.grid(axis='y', linestyle='--', alpha=0.7)


def _annotate_points(ax, scores, onset_ages, aliases, cluster_radius_px=12):
    """Label each point with its case alias, staggering labels for points that
    are at or near the same location (in pixel space) so they don't overlap."""
    points_px = ax.transData.transform(np.column_stack([scores, onset_ages]))

    # Assign each point to a coarse pixel-grid cell so points that are
    # identical or merely close together are clustered together.
    cell_of = [tuple(np.round(p / cluster_radius_px).astype(int)) for p in points_px]
    seen_in_cell = {}

    for (x, y), cell, alias in zip(zip(scores, onset_ages), cell_of, aliases):
        rank = seen_in_cell.get(cell, 0)
        seen_in_cell[cell] = rank + 1

        # Stagger successive labels in the same cluster further along a
        # diagonal so they fan out instead of stacking on top of each other.
        offset = (5 + rank * 10, 5 - rank * 6)
        if " " in alias:
            first_name, last_name = alias.split()
            alias_short = f"{first_name[0]}. {last_name}"
        else:
            alias_short = alias
        ax.annotate(alias_short, (x, y), fontsize=12, rotation=30 - 10 * rank,
                    xytext=offset, textcoords='offset points')


def _plot_scatter(ax: Axes, scores, onset_ages, haplotype_tested=None, aliases=None):
    """Scatter plot of onset ages vs scores, optionally colored by haplotype_tested
    and annotated with case aliases."""

    if haplotype_tested is not None:
        haplotype_tested = np.asarray(haplotype_tested)
        tested_mask = (haplotype_tested == 'yes')
        ax.scatter(scores[tested_mask], onset_ages[tested_mask],
                   alpha=0.7, color='tab:blue', label='haplotype tested')
        ax.scatter(scores[~tested_mask], onset_ages[~tested_mask],
                   alpha=0.7, color='tab:red', label='haplotype not tested')
        ax.legend()
    else:
        ax.scatter(scores, onset_ages, alpha=0.7)

    if aliases is not None:
        _annotate_points(ax, scores, onset_ages, aliases)

    ax.set_xlabel('Score', fontsize=24)
    ax.set_ylabel('Onset Age', fontsize=24)
    # ax.set_title('Onset Ages vs Scores')
    ax.tick_params(labelsize=18)
    ax.grid(True, alpha=0.3)


def _plot_violins_with_smart_data_clustering(ax: Axes, x: np.ndarray, y: np.ndarray):
    kde = estimate_density(x=x, bandwidth=0.2, prominence=0.03)
    labels = np.digitize(x, kde.boundaries)

    # Violins
    #
    positions, violins = [], []
    for cluster in np.unique(labels):
        mask = labels == cluster
        positions.append(np.median(x[mask]))
        violins.append(y[mask])

    vp = ax.violinplot(violins, positions=positions, widths=0.06, showmeans=False, showmedians=True, showextrema=False, )
    for body in vp["bodies"]: body.set_alpha(0.35)

    #
    # Cluster boundaries
    #
    for boundary in kde.boundaries:
        ax.axvline(boundary, ls="--", lw=1, color="gray", )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    # ax.set_title("Modal clustering using KDE minima")


def plot(scores, onset_ages, plot_type: str | None = None, haplotype_tested=None, aliases=None,
         save_to: str | None = None):
    """
    Plot onset ages vs scores.

    If plot_type is "whisker" or "violin": Creates a boxplot/violin plot where onset ages
                    are grouped by score bins (bins of width 0.1 from 0 to 1).
    If plot_type is None: Simple scatter plot. If haplotype_tested is provided, cases where
                    haplotype_tested != 'yes' are plotted in a different color. If aliases
                    is provided, each point is annotated with its case alias.
    If save_to is given, the figure is written to that file instead of being shown.
    """
    scores = np.asarray(scores)
    onset_ages = np.asarray(onset_ages)

    fig, ax = plt.subplots(figsize=(16, 12))
    if plot_type in ("whisker", "violin-basic"):
        _plot_boxplot_or_violin(ax, scores, onset_ages, plot_type)
    elif plot_type == "violin":
        _plot_violins_with_smart_data_clustering(ax, scores, onset_ages)
        _plot_scatter(ax, scores, onset_ages, haplotype_tested=haplotype_tested, aliases=None)
    else:
        _plot_scatter(ax, scores, onset_ages, haplotype_tested=haplotype_tested, aliases=aliases)

    if save_to:
        plt.savefig(save_to, bbox_inches='tight')
    else:
        plt.show()


def collect_onset_age_vs_score(cursor, publication_id: str, score_column: str):
    """Collect one point per case of the publication: the genotype score and the onset
    age, applying the same filters as the original 33_onset_age_vs_score.py - skip
    genotypes with a "bad" allele (a single variant with many gnomAD homozygotes) and
    homozygotes whose haplotype was not actually tested.

    Returns four parallel lists: scores, onset_ages, haplotype_tested, aliases."""
    qry = (f"select g.id, g.{score_column}, c.onset_age, c.haplotype_tested , c.patient_xref_id from genotypes g "
           "left join cases c on g.id=c.genotype_id "
           f"where  {score_column} is not null and onset_age > 0 "
           # "and haplotype_tested = 'yes' "
           f"and publication_id in ({publication_id}) "
           f"order by {score_column}")

    scores, onset_ages, haplotype_tested, aliases = [], [], [], []
    for genotype_id, score, onset_age, hapl_tested, patient_xref_id in error_intolerant_search(cursor, qry) or []:
        # bad allele has a single variant which has more than 10 homozygotes in gnomad
        if genotype_has_bad_allele(cursor, genotype_id): continue

        # filter suspicious
        qry = f"select allele_id1, allele_id2 from genotypes where id={genotype_id}"
        allele_id1, allele_id2 = hard_landing_search(cursor, qry)[0]
        homozygote = (allele_id1 == allele_id2)
        if homozygote and hapl_tested != 'yes':  continue

        scores.append(float(score))
        onset_ages.append(float(onset_age))
        haplotype_tested.append(hapl_tested)
        aliases.append(patient_xref_id)

    return scores, onset_ages, haplotype_tested, aliases


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=str(DEFAULT_DB),
                        help="path to the sqlite3 database (default: %(default)s)")
    parser.add_argument("-p", "--publication_id", type=str, default='70')
    parser.add_argument("-s", "--scoring_type", default="plain",
                        choices=['plain', 'score_w_dosage_compensation'])
    parser.add_argument("-v", "--visualization_type", default="scatter",
                        choices=['scatter', 'violin', 'violin-basic', 'whisker'])
    parser.add_argument("-o", "--out", default=None,
                        help="write the plot to this file instead of showing it interactively")
    return parser.parse_args()


def main():
    args = parse_args()
    score_column = "score" if args.scoring_type == "plain" else "score_w_dosage_compensation"

    db, cursor = connect(args.db)
    scores, onset_ages, haplotype_tested, aliases = \
        collect_onset_age_vs_score(cursor, args.publication_id, score_column)
    cursor.close()
    db.close()

    if not scores:
        print(f"no scored cases with onset age found for publication {args.publication_id} - "
              f"did you run score_n_store.py on {args.db} first?")
        return

    spearman = spearmanr(scores, onset_ages)
    pearson = pearsonr(scores, onset_ages)
    print(f"  Spearman correlation {spearman.statistic:.2f}   p-val: {spearman.pvalue:.2e}")
    print(f"  Pearson  correlation {pearson.statistic:.2f}    p-val: {pearson.pvalue:.2e}")

    aliases = []  # if this is empty the points are not annotated with aliases
    plot(scores, onset_ages, plot_type=args.visualization_type, haplotype_tested=haplotype_tested,
         aliases=aliases, save_to=args.out)


#########################
if __name__ == "__main__":
    main()
