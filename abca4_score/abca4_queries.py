"""SQLite counterpart of the original project's utils/abca4_mysql.py -
only the functions used by the functionality-score pipeline are ported.
The SQL and the logic are unchanged except that get_column_names lost its
database-name argument (meaningless in SQLite).
"""

from abca4_score.sqlite_utils import hard_landing_search, error_intolerant_search, get_column_names


def variant_ids_from_allele_id(cursor, allele_id) -> list[int]:
    qry = f"select variant_id from allele_variants where allele_id = {allele_id}"
    return [row[0] for row in hard_landing_search(cursor, qry)]


def allele_ids_containing_variant(cursor, variant_id) -> list[int]:
    # find all alleles that include variant_id among their variants (possibly along with others)
    qry = f"select distinct allele_id from allele_variants where variant_id = {variant_id}"
    return [row[0] for row in error_intolerant_search(cursor, qry) or []]


def allele_id_from_variant_ids(cursor, variant_ids) -> int | None:
    # find the allele whose set of variants is exactly variant_ids (order independent)
    variant_id_set = set(variant_ids)
    id_list = ",".join(str(v) for v in variant_id_set)
    qry = (f"select allele_id from allele_variants where variant_id in ({id_list}) "
           f"group by allele_id having count(*) = {len(variant_id_set)}")
    for row in error_intolerant_search(cursor, qry) or []:
        candidate_allele_id = row[0]
        qry = f"select count(*) from allele_variants where allele_id={candidate_allele_id}"
        total_variants = hard_landing_search(cursor, qry)[0][0]
        if total_variants == len(variant_id_set):
            return candidate_allele_id
    return None


def allele_ids_from_genotype_id(cursor, genotype_id) -> tuple[int, int]:
    qry = f"select allele_id1, allele_id2 from genotypes where id={genotype_id}"
    return hard_landing_search(cursor, qry)[0]


def score_from_genotype_id(cursor, genotype_id) -> float | None:
    qry = f"select score from genotypes where id={genotype_id}"
    return hard_landing_search(cursor, qry)[0][0]


def is_bad_allele(cursor, allele_id: int) -> bool:
    # bad = a single variant with more than 10 homozygotes in gnomAD
    variant_ids = variant_ids_from_allele_id(cursor, allele_id)
    varstr = ",".join(map(str, variant_ids))
    qry = f"select gnomad_homozygotes from variants where id in ({varstr})"
    homozygotes = hard_landing_search(cursor, qry)[0]
    return homozygotes is not None and all(h is not None and h > 10 for h in homozygotes)


def genotype_has_bad_allele(cursor, genotype_id: int) -> bool:
    allele_id1, allele_id2 = allele_ids_from_genotype_id(cursor, genotype_id)
    return is_bad_allele(cursor, allele_id1) or is_bad_allele(cursor, allele_id2)


def get_exp_characterization(cursor, allele_id, exp_header) -> list[dict]:
    qry = f"select * from exp_characterization where allele_id={allele_id}"
    exp = []
    for row in error_intolerant_search(cursor, qry) or []:
        exp.append(dict(zip(exp_header, row)))
    return exp


def related_allele_w_no_mrna(cursor, variant_ids: list[int]) -> int | None:
    exp_header = get_column_names(cursor, "exp_characterization")
    for variant_id in variant_ids:
        for other_allele_id in allele_ids_containing_variant(cursor, variant_id):
            if variant_ids_from_allele_id(cursor, other_allele_id) != [variant_id]: continue
            if variant_ids == [variant_id]: continue  # that is the same allele we started from
            exp_characterization = get_exp_characterization(cursor, other_allele_id, exp_header)
            if any(row['mRNA_level'] == 0 for row in exp_characterization):
                return other_allele_id
    return None


def allele_produces_no_mrna(cursor, variant_ids: list[int]) -> bool:
    # True if any variant in this allele also occurs, alone, in some other allele
    # that has been characterized and shows mRNA_level == 0 there
    return related_allele_w_no_mrna(cursor, variant_ids) is not None


def allele_has_early_stop(cursor, variant_ids: list[int]) -> bool:
    var_ids_string = ", ".join(str(v) for v in variant_ids)
    qry = f"select protein from variants where id in ({var_ids_string})"
    return any((row[0] is not None and 'Ter' in row[0]) for row in hard_landing_search(cursor, qry))


def get_gnomad_freqs(cursor, variant_ids: list[int]) -> dict[int, int]:
    qry = f"select id, gnomad_homozygotes from variants where id in ({','.join(str(v) for v in variant_ids)})"
    ret = error_intolerant_search(cursor, qry)
    gnomad_dict = dict(tuple(row) for row in ret) if ret else {}
    for vid in variant_ids:
        # -1 means: not found in gnomad at all (0 means that there are heterozygotes)
        if not gnomad_dict[vid]: gnomad_dict[vid] = -1
    return gnomad_dict


def find_corresponding_allele_wo_homozygotes(cursor, variant_ids: list[int]) -> int | None:
    gnomad_freq = get_gnomad_freqs(cursor, variant_ids)
    gnomad_homozygotes = [v for v in variant_ids if gnomad_freq[v] >= 5]
    if len(gnomad_homozygotes) < 1: return None
    if len(gnomad_homozygotes) == len(variant_ids): return None  # we should get rid of that some place else
    variant_not_homs = sorted(set(variant_ids).difference(set(gnomad_homozygotes)))

    return allele_id_from_variant_ids(cursor, variant_not_homs)
