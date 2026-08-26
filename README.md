# abca4_functionality_score

Standalone extraction of the ABCA4 genotype functionality-score pipeline, running on **sqlite3**.
The scoring pipeline itself uses the Python standard library only; the
onset-age-vs-score analysis additionally needs numpy, scipy and matplotlib.

For every genotype (a pair of alleles), the pipeline parametrizes each allele from
published experimental characterization data (expression, mRNA level, solubility,
ATPase activity, N-Ret-PE binding, ATR binding to ECD2 - each as a fraction of the
wild type), records how each allele's parametrization was arrived at, and combines
the two alleles into a 0..1 transport-competence score for the genotype. The
onset-age-vs-score analysis then correlates those scores with the reported age of
disease onset.

## Layout

```
score_n_store.py            entry point (port of 03_score_n_store.py)
onset_age_vs_score.py       onset age vs score correlation and plots
                            (port of abca4_60_production/generic/33_onset_age_vs_score.py)
abca4_score/
    sqlite_utils.py         sqlite port of the used parts of utils/mysql.py
    abca4_queries.py        sqlite port of the used parts of utils/abca4_mysql.py
    func_score_utils.py     port of func_score_utlis.py (the scoring logic, unchanged)
data/
    abca4_pub70_test.db     small sqlite3 test database (see below)
tests/
    test_scoring.py         end-to-end test of the scoring pipeline
    test_onset_age_vs_score.py  end-to-end test of the onset-age-vs-score analysis
    reference_scores.json   expected scores, computed by the ORIGINAL MySQL code
    reference_onset.json    expected onset-vs-score points and correlations, ditto
tools/
    extract_test_db.py      rebuilds data/ and tests/reference_*.json from MySQL
    (also the reference schema for a from-scratch database - see "Adding new data" below)
```

## Install

Requires Python 3.9+ (developed against 3.14) and `sqlite3`, which ships with
the standard library - `score_n_store.py` and the tests have no other
dependencies. `onset_age_vs_score.py` additionally needs numpy, scipy and
matplotlib:

```
git clone <this repo>
cd abca4_functionality_score
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt   # numpy, scipy, matplotlib - only needed for onset_age_vs_score.py
```

There is nothing to build or install into `site-packages`; the scripts are run
directly from the repo root as shown below.

## Usage

### 1. Quickstart

Score genotypes, originally from [Cobos et al](https://iovs.arvojournals.org/article.aspx?articleid=2811455),
in the provided database (in `data` directory), and plot the age of onset as a function of the scores
```
python3 score_n_store.py
python3 onset_age_vs_score.py
```
Note that both scripts proved usage statement when give the `--help` argument on the command line
```
python3 score_n_store.py --help
python3 onset_age_vs_score.py --help
```

### 2. Score the genotypes in a database

```
python3 score_n_store.py --db path/to/database.db [--assume_dosage_compensation] [--quiet]
```

Scores are written into `genotypes.score` (or `genotypes.score_w_dosage_compensation`
with the flag), and per-allele provenance into the `allele_scoring_source` table
(created on first run). The run modifies the database in place - work on a copy of
`data/abca4_pub70_test.db` if you want to keep it pristine. Run it once with each
flag combination if you want both score columns populated:

```
cp data/abca4_pub70_test.db  path/to/new/database/scored.db
python3 score_n_store.py --db path/to/new/database/scored.db
python3 score_n_store.py --db path/to/new/database/scored.db --assume_dosage_compensation
```

### 3. Plot age of onset vs score

Once a database has been scored (step 1), correlate the scores with the age of
onset recorded in `cases.onset_age`:

```
python3 onset_age_vs_score.py --db /tmp/scored.db [-p 70] \
        [-s score|score_w_dosage_compensation] \
        [-v scatter|violin|violin-basic|whisker] [-o plot.png]
```

- `--db` - path to the scored sqlite3 database (default: `data/abca4_pub70_test.db`)
- `-p/--publication-id` - restrict to one or more publication ids (default: `70`,
  matching the committed test db)
- `-s/--score-column` - which score column to plot against (default: `score`)
- `-v/--plot-type` - `scatter` (default), `violin` (KDE-clustered bins),
  `violin-basic`, or `whisker`
- `-o/--out` - write the figure to a file instead of opening an interactive window

This prints the Spearman and Pearson correlations and shows the plot (or writes
it to a file with `-o`, e.g. for use over ssh or in a script). As in the original,
genotypes with a "bad" allele (a single variant with more than 10 gnomAD
homozygotes) and homozygotes whose haplotype was not actually tested are excluded
from both the correlation and the plot.

## Adding new data to the sqlite3 database

The database has five tables (schema in `tools/extract_test_db.py:SQLITE_SCHEMA`):

| table                 | columns |
|------------------------|---|
| `variants`             | `id, cds, protein, gnomad_homozygotes` |
| `allele_variants`      | `allele_id, variant_id` (junction table: which variants make up an allele) |
| `exp_characterization` | `id, allele_id, publication_id, localization, mRNA_level, solubility, "N-Ret-PE_binding_no_ATP", "N-Ret-PE_binding_w_ATP", ATPase_activity_basal, ATPase_activity_retinal_stimulated, ATR_binding_to_ECD2, transport_between_liposomes` |
| `genotypes`            | `id, allele_id1, allele_id2, score, score_w_dosage_compensation` |
| `cases`                | `id, genotype_id, publication_id, patient_xref_id, onset_age, haplotype_tested` |

`allele_id`/`variant_id`/`genotype_id` are free-standing integers you choose
yourself (there is no separate `alleles` table - an allele *is* the set of
`allele_variants` rows sharing an `allele_id`); pick ids that don't already
exist in the target database. `exp_characterization.allele_id` is `NULL` for
the wild-type reference row of a publication. Score columns should be left
`NULL` on insert - they are filled in by `score_n_store.py`.

To add a new genotype end to end:

1. **Variants** - insert any variants that aren't already in the db, with their
   cDNA/protein change and gnomAD homozygote count (`gnomad_homozygotes` drives
   the "bad allele" filters, so leave it `NULL`/0 if unknown rather than guessing).
2. **Alleles** - for each allele, insert one `allele_variants` row per variant it
   carries (a simple missense allele is one row; a complex allele is several rows
   sharing the same `allele_id`).
3. **Experimental characterization** - insert `exp_characterization` rows for any
   allele (or the publication's wild-type row, `allele_id = NULL`) you have
   functional data for; leave a column `NULL` if that assay wasn't done. Alleles
   with no characterization at all are still handled by the pipeline's fallbacks
   (matching a related allele, stripping gnomAD-homozygote variants, etc. - see
   `abca4_score/func_score_utils.py`), just less precisely.
4. **Genotype** - insert one `genotypes` row per `(allele_id1, allele_id2)` pair,
   with the score columns `NULL`.
5. **Case** (optional, only needed for the onset-age analysis) - insert a `cases`
   row referencing the `genotype_id`, with `onset_age` and `haplotype_tested`
   (`'yes'`/`'no'`, relevant only for homozygous genotypes) filled in.

This can be done with the `sqlite3` CLI, or in Python with the same helpers the
pipeline uses:

```python
from abca4_score.sqlite_utils import connect

db, cursor = connect("path/to/database.db")
cursor.execute("insert into variants values (?,?,?,?)", (99001, "c.1000A>G", "p.Asn334Ser", 0))
cursor.execute("insert into allele_variants values (?,?)", (88001, 99001))
cursor.execute(
    "insert into exp_characterization "
    '(id, allele_id, publication_id, "mRNA_level", solubility) values (?,?,?,?,?)',
    (77001, 88001, 70, 0.8, 0.6))
cursor.execute("insert into genotypes (id, allele_id1, allele_id2) values (?,?,?)",
              (66001, 88001, 88001))
cursor.execute("insert into cases values (?,?,?,?,?,?)",
              (55001, 66001, 70, "new-patient-1", 12, "yes"))
db.close()
```

Then run `score_n_store.py` (twice, with and without `--assume_dosage_compensation`,
if you want both score columns) to score the new genotype.

To start a **brand-new, empty** database instead of adding to an existing one,
execute `tools.extract_test_db.SQLITE_SCHEMA` against a fresh sqlite3 file:

```python
import sqlite3
from tools.extract_test_db import SQLITE_SCHEMA

db = sqlite3.connect("path/to/new.db")
db.executescript(SQLITE_SCHEMA)
db.close()
```

## Tests

```
python3 -m unittest discover tests
```

`test_scoring.py` copies `data/abca4_pub70_test.db` to a temp file, runs
`score_n_store.py` on it in both dosage-compensation modes, and compares every
genotype score and every allele's scoring-source JSON against
`tests/reference_scores.json`, which was produced by running the *original*
MySQL-based code on the original database (at extraction time the reference also
matched the scores stored in MySQL for all 191 genotypes).

`test_onset_age_vs_score.py` scores a temp copy the same way, then checks that
the onset-age-vs-score point set (after the bad-allele and homozygote filters)
and the Spearman/Pearson correlations match `tests/reference_onset.json`, also
produced with the original MySQL-based code; finally it runs
`onset_age_vs_score.py` end to end and checks that it writes a plot.

Both test files can also be run individually, e.g. `python3 -m unittest
tests.test_onset_age_vs_score -v`, or with `pytest tests/` if pytest is installed.
