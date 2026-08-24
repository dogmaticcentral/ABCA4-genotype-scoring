#!/usr/bin/env python3
"""Port of the original abca4_90_functionality_score/03_score_n_store.py to SQLite.

For every genotype in the `genotypes` table, parametrize both alleles from the
experimental characterization data, record the provenance of each allele's
parametrization in `allele_scoring_source`, and, where a score can be computed,
store it in genotypes.score (or genotypes.score_w_dosage_compensation when
--assume_dosage_compensation is given).
"""
import argparse
from datetime import date
from pathlib import Path
from typing import Dict, Optional

from abca4_score.func_score_utils import AlleleScoringSource, NormalizedParametrization
from abca4_score.func_score_utils import competence_estimate, parametrize_genotype

from abca4_score.sqlite_utils import error_intolerant_search, connect, hard_landing_search, store_or_update

DEFAULT_DB = Path(__file__).parent / "data" / "abca4_pub70_test.db"


# we will take as a given that early stop codons get cleared through NMD
# TODO check for the cases that might evade NMD
def is_allele_ter(cursor, allele_id):
    qry = f"select protein from variants right join allele_variants on variants.id = variant_id where allele_id={allele_id}"
    proteins = [row[0] for row in hard_landing_search(cursor, qry) if row[0] is not None]
    return any('ter' in p.lower() for p in proteins)


def make_scoring_source_table_if_nonexistent(cursor):
    # the json object records how the allele was parametrized; scored_on is a date, no time needed
    qry = """
    CREATE TABLE IF NOT EXISTS allele_scoring_source (
        id         INTEGER PRIMARY KEY AUTOINCREMENT,
        allele_id  INTEGER NOT NULL UNIQUE,
        source     TEXT,
        scored_on  TEXT
    )
    """
    cursor.execute(qry)


def store_allele_scoring_source(cursor, allele_id: int, source: AlleleScoringSource):
    store_or_update(cursor, "allele_scoring_source",
                    fixed_fields={"allele_id": allele_id},
                    update_fields={"source": source.to_json(), "scored_on": date.today().isoformat()})


# a known problem. Gly1961Glu 10 homozygotes in gnomad
#     # in https://iovs.arvojournals.org/article.aspx?articleid=2680974  Molday says
#     # " G1961E mutation appears mild when this variant is retained in the membrane based on expression
#     #  studies and consistent with the relatively mild phenotype of individuals homozygous for this mutation.
#     #  However, after detergent solubilization, the variant is devoid of functional activity,
#     # including N-Ret-PE binding and ATPase activity. Detergent solubilization may adversely affect
#     # the functional activities of this mutant.
#     # I cannot find independent exp evaluation of this variant - let's not deal with that

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=str(DEFAULT_DB),
                        help="path to the sqlite3 database (default: %(default)s)")
    parser.add_argument("--assume_dosage_compensation", action="store_true", default=False,
                        help="assume dosage compensation when estimating competence")
    parser.add_argument("--quiet", action="store_true", default=False,
                        help="suppress the per-genotype printout")
    return parser.parse_args()


def main():
    args = parse_args()

    db, cursor = connect(args.db)
    make_scoring_source_table_if_nonexistent(cursor)

    qry = "select id, allele_id1, allele_id2 from genotypes "
    for genotype_id, allele_id1, allele_id2 in hard_landing_search(cursor, qry):
        genotype_params: Dict[int, Optional[NormalizedParametrization]]
        genotype_sources: Dict[int, AlleleScoringSource]
        genotype_params, genotype_sources = parametrize_genotype(cursor, (allele_id1, allele_id2), verbose=False)

        # record how each allele was (or was not) parametrized, regardless of whether we end up scoring
        for allele_id, source in genotype_sources.items():
            store_allele_scoring_source(cursor, allele_id, source)

        if any(prms is None for prms in genotype_params.values()):
            continue
        score = competence_estimate(genotype_params, is_homozygote=(allele_id1 == allele_id2),
                                    assume_dosage_compensation=args.assume_dosage_compensation,
                                    verbose=not args.quiet)
        if score < 0: continue  # no info to calculate the score - not sure if this can happen in this version of the code
        score_column = "score_w_dosage_compensation" if args.assume_dosage_compensation else "score"
        qry = f"update genotypes set {score_column}={score} where id={genotype_id}"
        error_intolerant_search(cursor, qry)

    cursor.close()
    db.close()


################################
if __name__ == "__main__":
    main()
