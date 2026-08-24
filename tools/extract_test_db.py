#!/usr/bin/env python3
"""Build the small SQLite test database (data/abca4_pub<N>_test.db) from the
original abca4 MySQL database, together with the reference scores
(tests/reference_scores.json) and the reference onset-age-vs-score dataset
(tests/reference_onset.json) against which the SQLite port is tested.

The extraction is read-only on the MySQL side. It needs the original abca4
repo (for its utils.mysql and the original scoring code) and a working
~/.abca4_conf, so it is meant to be run in the environment where the MySQL
database lives:

    python3 tools/extract_test_db.py [--abca4-repo /path/to/abca4] [--publication-id N]

What is copied - the minimal closure that the scoring pipeline and the
onset-age-vs-score analysis can touch for the cases of the chosen publication
(default 70):

  cases              only the chosen publication_id; the identifying columns
                     plus onset_age and haplotype_tested (needed by the
                     onset-age-vs-score analysis)
  genotypes          genotypes of those cases; score columns left NULL
                     (the committed db is the pristine, unscored state)
  allele_variants    rows for: the alleles of those genotypes ("primary"
                     alleles), every single-variant allele carrying one of the
                     primary alleles' variants (needed by the
                     related_allele_w_no_mrna / individual-variants fallbacks),
                     and every allele matching a primary allele stripped of its
                     gnomAD-homozygote variants (needed by
                     find_corresponding_allele_wo_homozygotes)
  variants           variants of the alleles above; only id, cds, protein,
                     gnomad_homozygotes (the columns the pipeline reads,
                     plus cds for human readability)
  exp_characterization  rows for the alleles above, plus the wild-type rows
                     (allele_id is null) of every publication those rows cite;
                     the free-text notes column is dropped

The reference scores are computed with the ORIGINAL MySQL-based code
(abca4_90_functionality_score.func_score_utlis) run against MySQL, for both
dosage-compensation settings, so the test asserts that the SQLite port working
off the extracted db reproduces the original pipeline exactly. Likewise the
reference onset dataset is built with the original filtering code
(utils.abca4_mysql.genotype_has_bad_allele) run against MySQL.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PUBLICATION_ID = 70
GNOMAD_HOMOZYGOTE_CUTOFF = 5  # same cutoff as find_corresponding_allele_wo_homozygotes

EXP_COLUMNS = ["id", "allele_id", "publication_id", "localization", "mRNA_level", "solubility",
               "N-Ret-PE_binding_no_ATP", "N-Ret-PE_binding_w_ATP", "ATPase_activity_basal",
               "ATPase_activity_retinal_stimulated", "ATR_binding_to_ECD2", "transport_between_liposomes"]

SQLITE_SCHEMA = """
CREATE TABLE cases (
    id                INTEGER PRIMARY KEY,
    genotype_id       INTEGER,
    publication_id    INTEGER,
    patient_xref_id   TEXT,
    onset_age         INTEGER,
    haplotype_tested  TEXT
);
CREATE TABLE genotypes (
    id          INTEGER PRIMARY KEY,
    allele_id1  INTEGER NOT NULL,
    allele_id2  INTEGER NOT NULL,
    score       REAL,
    score_w_dosage_compensation REAL
);
CREATE TABLE allele_variants (
    allele_id   INTEGER NOT NULL,
    variant_id  INTEGER NOT NULL,
    PRIMARY KEY (allele_id, variant_id)
);
CREATE TABLE variants (
    id                 INTEGER PRIMARY KEY,
    cds                TEXT,
    protein            TEXT,
    gnomad_homozygotes INTEGER
);
CREATE TABLE exp_characterization (
    id              INTEGER PRIMARY KEY,
    allele_id       INTEGER,
    publication_id  INTEGER NOT NULL,
    localization    TEXT,
    mRNA_level      REAL,
    solubility      REAL,
    "N-Ret-PE_binding_no_ATP"  REAL,
    "N-Ret-PE_binding_w_ATP"   REAL,
    ATPase_activity_basal      REAL,
    ATPase_activity_retinal_stimulated REAL,
    ATR_binding_to_ECD2        REAL,
    transport_between_liposomes REAL
);
"""


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abca4-repo", default="/home/ivana/academia/projects/abca4",
                        help="path to the original abca4 repo (default: %(default)s)")
    parser.add_argument("--publication-id", type=int, default=DEFAULT_PUBLICATION_ID,
                        help="publication whose cases the test db is built from (default: %(default)s)")
    parser.add_argument("--out-db", default=None,
                        help="default: data/abca4_pub<N>_test.db")
    parser.add_argument("--out-ref", default=str(REPO_ROOT / "tests" / "reference_scores.json"))
    parser.add_argument("--out-onset-ref", default=str(REPO_ROOT / "tests" / "reference_onset.json"))
    args = parser.parse_args()
    if args.out_db is None:
        args.out_db = str(REPO_ROOT / "data" / f"abca4_pub{args.publication_id}_test.db")
    return args


def variant_ids_of_allele(cursor, hard_landing_search, allele_id):
    qry = f"select variant_id from allele_variants where allele_id = {allele_id}"
    return [row[0] for row in hard_landing_search(cursor, qry)]


def collect_closure(cursor, hard_landing_search, error_intolerant_search, publication_id):
    # the cases and their genotypes
    qry = (f"select id, genotype_id, publication_id, patient_xref_id, onset_age, haplotype_tested "
           f"from cases where publication_id = {publication_id}")
    case_rows = hard_landing_search(cursor, qry)
    genotype_ids = sorted({row[1] for row in case_rows if row[1] is not None})

    qry = f"select id, allele_id1, allele_id2 from genotypes where id in ({','.join(map(str, genotype_ids))})"
    genotype_rows = hard_landing_search(cursor, qry)

    primary_alleles = sorted({a for _, a1, a2 in genotype_rows for a in (a1, a2)})

    # variants of the primary alleles
    primary_variants = set()
    variants_by_allele = {}
    for allele_id in primary_alleles:
        vids = variant_ids_of_allele(cursor, hard_landing_search, allele_id)
        variants_by_allele[allele_id] = vids
        primary_variants.update(vids)

    allele_ids = set(primary_alleles)

    # related single-variant alleles: for each primary variant, every allele whose
    # only variant it is (the related_allele_w_no_mrna / individual-variants fallbacks)
    vstr = ",".join(map(str, sorted(primary_variants)))
    qry = f"select distinct allele_id from allele_variants where variant_id in ({vstr})"
    for row in hard_landing_search(cursor, qry):
        candidate = row[0]
        qry = f"select count(*) from allele_variants where allele_id = {candidate}"
        if hard_landing_search(cursor, qry)[0][0] == 1:
            allele_ids.add(candidate)

    # alleles corresponding to a primary allele stripped of its gnomAD-homozygote
    # variants (find_corresponding_allele_wo_homozygotes)
    for allele_id in primary_alleles:
        vids = variants_by_allele[allele_id]
        qry = f"select id, gnomad_homozygotes from variants where id in ({','.join(map(str, vids))})"
        homozygotes = {row[0]: (row[1] if row[1] else -1) for row in hard_landing_search(cursor, qry)}
        hom_variants = [v for v in vids if homozygotes[v] >= GNOMAD_HOMOZYGOTE_CUTOFF]
        if not hom_variants or len(hom_variants) == len(vids):
            continue
        subset = sorted(set(vids).difference(hom_variants))
        sstr = ",".join(map(str, subset))
        qry = (f"select allele_id from allele_variants where variant_id in ({sstr}) "
               f"group by allele_id having count(*) = {len(subset)}")
        for row in error_intolerant_search(cursor, qry) or []:
            candidate = row[0]
            qry = f"select count(*) from allele_variants where allele_id = {candidate}"
            if hard_landing_search(cursor, qry)[0][0] == len(subset):
                allele_ids.add(candidate)

    # allele_variants rows and the full variant set for the collected alleles
    allele_variant_rows = []
    variant_ids = set()
    for allele_id in sorted(allele_ids):
        for vid in variant_ids_of_allele(cursor, hard_landing_search, allele_id):
            allele_variant_rows.append((allele_id, vid))
            variant_ids.add(vid)

    qry = (f"select id, cds, protein, gnomad_homozygotes from variants "
           f"where id in ({','.join(map(str, sorted(variant_ids)))})")
    variant_rows = hard_landing_search(cursor, qry)

    # experimental characterization of the collected alleles + the wild-type rows
    # of every publication involved
    col_str = ", ".join(f"`{c}`" for c in EXP_COLUMNS)
    astr = ",".join(map(str, sorted(allele_ids)))
    qry = f"select {col_str} from exp_characterization where allele_id in ({astr})"
    exp_rows = [tuple(row) for row in error_intolerant_search(cursor, qry) or []]
    exp_pub_ids = sorted({row[2] for row in exp_rows})
    if exp_pub_ids:
        qry = (f"select {col_str} from exp_characterization where allele_id is null "
               f"and publication_id in ({','.join(map(str, exp_pub_ids))})")
        exp_rows += [tuple(row) for row in hard_landing_search(cursor, qry)]

    return {
        "cases": [tuple(row) for row in case_rows],
        "genotypes": [tuple(row) for row in genotype_rows],
        "allele_variants": allele_variant_rows,
        "variants": [tuple(row) for row in variant_rows],
        "exp_characterization": exp_rows,
    }


def write_sqlite(out_db: Path, closure: dict):
    out_db.parent.mkdir(parents=True, exist_ok=True)
    if out_db.exists():
        out_db.unlink()
    db = sqlite3.connect(out_db)
    db.executescript(SQLITE_SCHEMA)
    db.executemany("insert into cases values (?,?,?,?,?,?)", closure["cases"])
    # score columns deliberately NULL - the committed db is the unscored state
    db.executemany("insert into genotypes (id, allele_id1, allele_id2) values (?,?,?)",
                   closure["genotypes"])
    db.executemany("insert into allele_variants values (?,?)", closure["allele_variants"])
    db.executemany("insert into variants values (?,?,?,?)", closure["variants"])
    placeholders = ",".join("?" for _ in EXP_COLUMNS)
    db.executemany(f"insert into exp_characterization values ({placeholders})",
                   closure["exp_characterization"])
    db.commit()
    db.close()


def compute_reference_scores(cursor, genotype_rows):
    """Run the ORIGINAL scoring code (against MySQL) for the extracted genotypes."""
    from abca4_90_functionality_score.func_score_utlis import competence_estimate, parametrize_genotype

    reference = {}
    for genotype_id, allele_id1, allele_id2 in genotype_rows:
        params, sources = parametrize_genotype(cursor, (allele_id1, allele_id2), verbose=False)
        entry = {"score": None, "score_w_dosage_compensation": None,
                 "scoring_source": {str(allele_id): json.loads(source.to_json())
                                    for allele_id, source in sources.items()}}
        if all(prms is not None for prms in params.values()):
            for key, dosage in (("score", False), ("score_w_dosage_compensation", True)):
                score = competence_estimate(params, is_homozygote=(allele_id1 == allele_id2),
                                            assume_dosage_compensation=dosage, verbose=False)
                if score >= 0:
                    entry[key] = score
        reference[str(genotype_id)] = entry
    return reference


def compute_reference_onset(cursor, hard_landing_search, reference, publication_id):
    """Build the reference for the onset-age-vs-score analysis: for both score columns,
    the (alias, onset_age, score) points that survive the original filtering of
    33_onset_age_vs_score.py - using the ORIGINAL genotype_has_bad_allele (run against
    MySQL) and the freshly computed reference scores - plus the resulting Spearman and
    Pearson correlations."""
    from scipy.stats import spearmanr, pearsonr
    from utils.abca4_mysql import genotype_has_bad_allele

    qry = (f"select genotype_id, onset_age, haplotype_tested, patient_xref_id from cases "
           f"where publication_id = {publication_id} and genotype_id is not null and onset_age > 0")
    case_rows = hard_landing_search(cursor, qry)

    onset_reference = {}
    for score_key in ("score", "score_w_dosage_compensation"):
        points = []
        for genotype_id, onset_age, hapl_tested, alias in case_rows:
            score = reference[str(genotype_id)][score_key]
            if score is None: continue
            if genotype_has_bad_allele(cursor, genotype_id): continue
            qry = f"select allele_id1, allele_id2 from genotypes where id = {genotype_id}"
            allele_id1, allele_id2 = hard_landing_search(cursor, qry)[0]
            if allele_id1 == allele_id2 and hapl_tested != 'yes': continue
            points.append([alias, float(onset_age), score])
        points.sort(key=lambda p: (p[0] or "", p[1], p[2]))
        scores = [p[2] for p in points]
        onset_ages = [p[1] for p in points]
        spearman = spearmanr(scores, onset_ages)
        pearson = pearsonr(scores, onset_ages)
        onset_reference[score_key] = {
            "points": points,
            "spearman": spearman.statistic, "spearman_pvalue": spearman.pvalue,
            "pearson": pearson.statistic, "pearson_pvalue": pearson.pvalue,
        }
        print(f"onset-age-vs-score reference ({score_key}): {len(points)} points, "
              f"spearman {spearman.statistic:.2f}, pearson {pearson.statistic:.2f}")
    return onset_reference


def sanity_check_against_stored_scores(cursor, hard_landing_search, reference):
    """Compare the freshly computed reference against the score column already
    stored in MySQL (decimal(5,2), so compare after rounding)."""
    mismatches = 0
    for genotype_id, entry in reference.items():
        qry = f"select score from genotypes where id = {genotype_id}"
        stored = hard_landing_search(cursor, qry)[0][0]
        computed = entry["score"]
        if stored is None and computed is None:
            continue
        if stored is None or computed is None or abs(float(stored) - round(computed, 2)) > 0.011:
            print(f"  genotype {genotype_id}: stored score {stored}, freshly computed "
                  f"{None if computed is None else round(computed, 2)}")
            mismatches += 1
    print(f"sanity check vs scores stored in mysql: {mismatches} mismatch(es) "
          f"out of {len(reference)} genotypes")


def main():
    args = parse_args()
    sys.path.insert(0, args.abca4_repo)
    from utils.mysql import abca4_connect, hard_landing_search, error_intolerant_search

    db, cursor = abca4_connect()

    closure = collect_closure(cursor, hard_landing_search, error_intolerant_search, args.publication_id)
    for table, rows in closure.items():
        print(f"{table}: {len(rows)} rows")

    write_sqlite(Path(args.out_db), closure)
    print(f"wrote {args.out_db}")

    reference = compute_reference_scores(cursor, closure["genotypes"])
    out_ref = Path(args.out_ref)
    out_ref.parent.mkdir(parents=True, exist_ok=True)
    with open(out_ref, "w") as f:
        json.dump(reference, f, indent=2)
    print(f"wrote {out_ref}")

    onset_reference = compute_reference_onset(cursor, hard_landing_search, reference, args.publication_id)
    with open(args.out_onset_ref, "w") as f:
        json.dump(onset_reference, f, indent=2)
    print(f"wrote {args.out_onset_ref}")

    sanity_check_against_stored_scores(cursor, hard_landing_search, reference)

    cursor.close()
    db.close()


if __name__ == "__main__":
    main()
