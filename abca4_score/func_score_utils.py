"""Port of the original abca4_90_functionality_score/func_score_utlis.py to SQLite.

The scoring logic is untouched; the only changes are the imports (the SQLite
utility modules instead of the MySQL ones), the adjusted get_column_names
signature, and the removal of a few functions that the score-and-store
pipeline never calls (homozygotes_exist, given, get_raw_* helpers).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, fields
from typing import Optional, List, Dict, Tuple

from abca4_score.abca4_queries import find_corresponding_allele_wo_homozygotes, allele_has_early_stop, \
    variant_ids_from_allele_id, related_allele_w_no_mrna, allele_id_from_variant_ids
from abca4_score.sqlite_utils import get_column_names, error_intolerant_search, hard_landing_search

quantitative_attrs = ['ATPase_activity', 'ATR_binding_to_ECD2', 'N-Ret-PE_binding', 'mRNA_level', 'expression']


@dataclass
class RawParametrization:
    """One row of the `exp_characterization` table (measurements stored as pct of the wild type).

    Attribute names mirror the DB columns, except that the two hyphenated columns
    `N-Ret-PE_binding_no_ATP` / `N-Ret-PE_binding_w_ATP` become
    `N_Ret_PE_binding_no_ATP` / `N_Ret_PE_binding_w_ATP`, since hyphens are not
    valid in Python identifiers. The `notes` column is intentionally dropped.
    """
    publication_id: int
    id: Optional[int] = None
    allele_id: Optional[int] = None
    localization: Optional[str] = None
    mRNA_level: Optional[float] = None
    solubility: Optional[float] = None
    N_Ret_PE_binding_no_ATP: Optional[float] = None
    N_Ret_PE_binding_w_ATP: Optional[float] = None
    ATPase_activity_basal: Optional[float] = None
    ATPase_activity_retinal_stimulated: Optional[float] = None
    ATR_binding_to_ECD2: Optional[float] = None
    transport_between_liposomes: Optional[float] = None
    corresponding_wt_params: Optional[RawParametrization] = None

    @classmethod
    def from_dict(cls, row: dict) -> RawParametrization:
        """Build from a `column_name -> value` dict, mapping hyphens to underscores
        and ignoring columns (e.g. `notes`) that are not fields of this class."""
        field_names = {f.name for f in fields(cls)}
        return cls(**{k.replace("-", "_"): v for k, v in row.items() if k.replace("-", "_") in field_names})


@dataclass
class NormalizedParametrization:
    """Functional parameters expressed as a fraction of the wild type (clamped to 0..1).

    Field names match the keys of the historical `missense_parametrization` dict,
    with `N-Ret-PE_binding` spelled `N_Ret_PE_binding`. Each `*_is_assumed` flag
    records whether the accompanying value was actually measured (False) or guessed
    (True). Defaults describe an assumed, wild-type-like allele.
    """
    ATPase_activity: float = 1.0
    ATPase_activity_is_assumed: bool = True
    N_Ret_PE_binding: float = 1.0
    N_Ret_PE_binding_is_assumed: bool = True
    ATR_binding_to_ECD2: float = 1.0
    ATR_binding_to_ECD2_is_assumed: bool = True
    mRNA_level: float = 1.0
    mRNA_level_is_assumed: bool = True
    solubility: float = 1.0
    solubility_is_assumed: bool = True
    expression: float = 1.0
    expression_is_assumed: bool = True


# a missense change is assumed intact everywhere except that we take expression at face value (0, measured)
missense_parametrization = NormalizedParametrization(
    solubility_is_assumed=False,
    expression=0.0,
    expression_is_assumed=False,
)


def early_stop_parametrization() -> NormalizedParametrization:
    # early stop: no full-length mRNA, hence no protein expressed
    return NormalizedParametrization(
        mRNA_level=0.0, mRNA_level_is_assumed=False,
        expression=0.0, expression_is_assumed=False,
    )


@dataclass
class AlleleScoringSource:
    """Records how an allele's functional parametrization was arrived at, so the decision
    can be reconstructed later. Serialized to JSON and stored in `allele_scoring_source`.

    `method` is one of:
      "none"                                    - nothing could be deduced (serializes to {})
      "early_stop"                              - the allele's protein carries a premature stop
      "direct"                                  - the allele itself is experimentally characterized
      "related_allele_no_mrna"                  - a variant of this allele, alone in another allele, shows mRNA == 0
      "characterized_allele_wo_common_variants" - common (homozygous) variants stripped, remainder is characterized
      "individual_variants"                     - built from the per-variant minimum of characterized single variants

    The extra fields are populated only when relevant:
      `exp_characterization_ids` - ids in `exp_characterization` the estimate rests on
      `allele_id`                - the other allele whose characterization was borrowed
      `protein`                  - the Ter-carrying protein descriptor(s) (for "early_stop")
      `variants`                 - per-variant trace (for "individual_variants"), each entry a dict with
                                   `variant_id` and either `early_stop: True` or
                                   `characterized_in_allele_id` + `exp_characterization_ids`
    """
    method: str = "none"
    exp_characterization_ids: Optional[List[int]] = None
    allele_id: Optional[int] = None
    protein: Optional[List[str]] = None
    variants: Optional[List[dict]] = None

    def to_json(self) -> str:
        if self.method == "none":
            return json.dumps({})
        payload: Dict = {"method": self.method}
        for key in ("exp_characterization_ids", "allele_id", "protein", "variants"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = value
        return json.dumps(payload)


def normalize(prm: RawParametrization) -> NormalizedParametrization:
    normalized = NormalizedParametrization()
    wt = prm.corresponding_wt_params

    # normalized expression: E-factor in  https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC7755071/
    for attribute in ('mRNA_level', 'solubility', 'ATR_binding_to_ECD2'):
        prm_val = getattr(prm, attribute)
        if prm_val is None:
            setattr(normalized, attribute, 1)
            setattr(normalized, attribute + "_is_assumed", True)
        else:
            setattr(normalized, attribute, prm_val / getattr(wt, attribute))
            setattr(normalized, attribute + "_is_assumed", False)

    if normalized.mRNA_level_is_assumed == normalized.solubility_is_assumed:
        # normalized.expression = min(normalized.mRNA_level, normalized.solubility)
        normalized.expression = normalized.mRNA_level * normalized.solubility
        normalized.expression_is_assumed = normalized.mRNA_level_is_assumed

    elif not normalized.mRNA_level_is_assumed:
        normalized.expression = normalized.mRNA_level
        normalized.expression_is_assumed = False

    else:
        normalized.expression = normalized.solubility
        normalized.expression_is_assumed = False

    # the more complicated ones
    if prm.ATPase_activity_basal is None and prm.ATPase_activity_retinal_stimulated is None:
        normalized.ATPase_activity = 1
        normalized.ATPase_activity_is_assumed = True


    elif prm.ATPase_activity_basal is not None and prm.ATPase_activity_retinal_stimulated is not None:
        diff_mutant = prm.ATPase_activity_retinal_stimulated - prm.ATPase_activity_basal
        diff_wt = wt.ATPase_activity_retinal_stimulated - wt.ATPase_activity_basal
        # diff_mutant = prm.ATPase_activity_basal
        # diff_wt =  wt.ATPase_activity_basal
        normalized.ATPase_activity = diff_mutant / diff_wt

    # ATPase assay just measures the release of the phosphate:  https://en.wikipedia.org/wiki/ATPase_assay
    # "ATP hydrolysis yields inorganic phosphate (Pi), which can be measured by a simple colorimetric reaction.
    # The amount of Pi liberated is directly proportional to the activity of the transporter."
    # so even though the retinal is not making much difference. the transporter is not as dead as it could be;
    # see for example here https://www.nature.com/articles/ng1000_242#Sec2, Fig 4 in particular
    # T1526M in particular - it is around 50% of the WT and stays so even in the presence of retinal, as
    # documented also in https://open.library.ubc.ca/media/stream/pdf/831/1.0090157/1 , page 90 (there actually
    # the activity is close to WT when unstimulated)
    # however, referring again to https://www.nature.com/articles/ng1000_242#Sec2,
    # figure 5, tha ATPase activity can be flat out 0
    elif prm.ATPase_activity_retinal_stimulated is not None:
        normalized.ATPase_activity = prm.ATPase_activity_retinal_stimulated / wt.ATPase_activity_retinal_stimulated
        normalized.ATPase_activity_is_assumed = False

    elif prm.ATPase_activity_basal is not None:
        normalized.ATPase_activity = prm.ATPase_activity_basal / wt.ATPase_activity_basal
        normalized.ATPase_activity_is_assumed = False

    else:
        raise RuntimeError("unreachable ATPase_activity branch")

    if prm.N_Ret_PE_binding_no_ATP is None and prm.N_Ret_PE_binding_w_ATP is None:
        normalized.N_Ret_PE_binding = 1
        normalized.N_Ret_PE_binding_is_assumed = True

    elif prm.N_Ret_PE_binding_no_ATP is not None and prm.N_Ret_PE_binding_w_ATP is None:
        normalized.N_Ret_PE_binding = prm.N_Ret_PE_binding_no_ATP / wt.N_Ret_PE_binding_no_ATP
        normalized.N_Ret_PE_binding_is_assumed = False

    elif prm.N_Ret_PE_binding_no_ATP is not None and prm.N_Ret_PE_binding_w_ATP is not None:
        # S-factor in https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC7755071/
        diff_mutant = prm.N_Ret_PE_binding_no_ATP - prm.N_Ret_PE_binding_w_ATP
        diff_wt = wt.N_Ret_PE_binding_no_ATP - wt.N_Ret_PE_binding_w_ATP
        normalized.N_Ret_PE_binding = diff_mutant / diff_wt
        normalized.N_Ret_PE_binding_is_assumed = False
    else:
        raise RuntimeError("unreachable N-Ret-PE_binding branch")

    for attribute in quantitative_attrs:
        attr = attribute.replace("-", "_")
        setattr(normalized, attr, min(1, max(0, getattr(normalized, attr))))

    return normalized


def get_params_for_characterized_allele(cursor, allele_id: int | List) -> List:
    exp_header = get_column_names(cursor, 'exp_characterization')
    if isinstance(allele_id, int):
        qry = f"select * from exp_characterization where allele_id = {allele_id} "
    else:
        allele_id_str = ", ".join(map(str, allele_id))
        qry = f"select * from exp_characterization where allele_id in ({allele_id_str}) "

    allele_params = []
    for line in hard_landing_search(cursor, qry):
        named_values = dict(zip(exp_header, line))
        raw_params = RawParametrization.from_dict(named_values)
        allele_params.append(raw_params)
        # the corresponding wt values should be retrievable from the database
        pub_id = named_values['publication_id']
        qry = f"select * from exp_characterization where allele_id is null and publication_id = {pub_id} "
        ret_wt = error_intolerant_search(cursor, qry) or []
        if len(ret_wt) == 0:
            raise Exception(f"No wt values found for publication id {pub_id}")
        elif len(ret_wt) > 1:
            raise Exception(f"Multiple wt values found for publication id {pub_id}")
        named_values = dict(zip(exp_header, ret_wt[0]))
        # TODO maybe this could be optimized a bit by filling wt for all wt
        # in a separate pass - in case they are from the same publication
        raw_params.corresponding_wt_params = RawParametrization.from_dict(named_values)
    return allele_params


def norm_allele_prms_avgd_over_sources(raw_allele_params: list[
    RawParametrization], verbose: bool) -> NormalizedParametrization | None:
    if len(raw_allele_params) == 0:
        return None

    norm_params = [normalize(prm) for prm in raw_allele_params]
    if verbose:
        print("\nparameters raw:")
        print(raw_allele_params)
        print("\nparameters normalized:")
        print(norm_params)
        print()

    avg_params = NormalizedParametrization()
    #  do not average if a value was assumed, rather than given in the publication
    for q in quantitative_attrs:
        attr = q.replace("-", "_")
        not_assumed = [getattr(prm, attr) for prm in norm_params if not getattr(prm, attr + "_is_assumed")]
        if len(not_assumed) > 0:
            setattr(avg_params, q, sum(not_assumed) / len(not_assumed))
        else:  # if everybody's guessing, then we can take the average just the same, what the heck
            all_guessed = [getattr(prm, attr) for prm in norm_params]
            setattr(avg_params, q, sum(all_guessed) / len(all_guessed))
    return avg_params


def allele_is_characterized(cursor, allele_id: int) -> bool:
    qry = f"select count(*) from exp_characterization where allele_id = {allele_id} "
    return hard_landing_search(cursor, qry)[0][0] > 0


def exp_characterization_ids(raw_params: List[RawParametrization]) -> List[int]:
    return [prm.id for prm in raw_params if prm.id is not None]


def ter_protein_tokens(cursor, variant_ids: List[int]) -> List[str]:
    if not variant_ids:
        return []
    var_ids_string = ", ".join(str(v) for v in variant_ids)
    qry = f"select protein from variants where id in ({var_ids_string}) "
    return [row[0] for row in hard_landing_search(cursor, qry) if row[0] is not None and 'Ter' in row[0]]


def min_over_normalized(params_list: List[NormalizedParametrization]) -> NormalizedParametrization:
    """Combine several normalized parametrizations into one by taking, for each characteristic,
    the minimum over the inputs. Used when an allele carrying several variants is estimated from
    the individual variants: the allele is assumed to be at least as impaired as its worst variant."""
    if len(params_list) == 0:
        raise Exception("min_over_normalized called with an empty list")

    combined = NormalizedParametrization()
    for attribute in quantitative_attrs:
        attr = attribute.replace("-", "_")
        setattr(combined, attr, min(getattr(prm, attr) for prm in params_list))
        # the characteristic counts as measured for the allele if it was measured for at least one variant
        all_assumed = all(getattr(prm, attr + "_is_assumed") for prm in params_list)
        setattr(combined, attr + "_is_assumed", all_assumed)
    return combined


def characterized_single_variant_allele(cursor, variant_id: int) -> Optional[int]:
    """Return the id of the allele whose only variant is `variant_id`, if that allele exists
    and has been experimentally characterized; otherwise None."""
    allele_id = allele_id_from_variant_ids(cursor, [variant_id])
    if allele_id is None:
        return None
    if not allele_is_characterized(cursor, allele_id):
        return None
    return allele_id


def parametrize_from_individual_variants(
        cursor, allele_id: int, variant_ids: List[int], verbose: bool
) -> Tuple[Optional[NormalizedParametrization], AlleleScoringSource]:
    """An allele with several variants and no characterization of its own
    is estimated from the individual variants. Each variant is characterized separately - as an early
    stop, or through an allele that contains that single variant - and the allele is assigned the
    minimum of each characteristic (see min_over_normalized)."""
    per_variant_params: List[NormalizedParametrization] = []
    trace: List[dict] = []
    for variant_id in variant_ids:
        if allele_has_early_stop(cursor, [variant_id]):
            per_variant_params.append(early_stop_parametrization())
            trace.append({"variant_id": variant_id, "early_stop": True})
            continue

        char_allele_id = characterized_single_variant_allele(cursor, variant_id)
        if char_allele_id is not None:
            raw = get_params_for_characterized_allele(cursor, char_allele_id)
            norm = norm_allele_prms_avgd_over_sources(raw, verbose)
            if norm is not None:
                per_variant_params.append(norm)
                trace.append({"variant_id": variant_id,
                              "characterized_in_allele_id": char_allele_id,
                              "exp_characterization_ids": exp_characterization_ids(raw)})
                continue

        # nothing found for this variant - it contributes no lower bound (i.e. is taken as wild-type)
        trace.append({"variant_id": variant_id})

    if len(per_variant_params) == 0:
        return None, AlleleScoringSource(method="none")

    return min_over_normalized(per_variant_params), AlleleScoringSource(method="individual_variants", variants=trace)


def parametrize_allele(cursor, allele_id: int, verbose: bool = False
                       ) -> Tuple[Optional[NormalizedParametrization], AlleleScoringSource]:
    """Normalized (and averaged over sources) parametrization of a single allele, together with a record
    of how it was obtained. Returns (None, AlleleScoringSource(method="none")) when nothing can be deduced."""
    variant_ids = variant_ids_from_allele_id(cursor, allele_id)

    # (1) the allele itself has been characterized in the literature
    if allele_is_characterized(cursor, allele_id):
        raw = get_params_for_characterized_allele(cursor, allele_id)
        source = AlleleScoringSource(method="direct", exp_characterization_ids=exp_characterization_ids(raw))
        return norm_allele_prms_avgd_over_sources(raw, verbose), source

    # (2) the protein carries a premature stop codon - no functional product
    if allele_has_early_stop(cursor, variant_ids):
        source = AlleleScoringSource(method="early_stop", protein=ter_protein_tokens(cursor, variant_ids))
        return early_stop_parametrization(), source

    # (3) fall back on the individual variants:
    #     (a) a variant of this allele, on its own in another allele, abolishes the mRNA
    related_no_mrna_allele_id = related_allele_w_no_mrna(cursor, variant_ids)
    if related_no_mrna_allele_id is not None:
        raw = get_params_for_characterized_allele(cursor, related_no_mrna_allele_id)
        source = AlleleScoringSource(method="related_allele_no_mrna", allele_id=related_no_mrna_allele_id,
                                     exp_characterization_ids=exp_characterization_ids(raw))
        return norm_allele_prms_avgd_over_sources(raw, verbose), source

    #     (b) strip the common (homozygous) variants and use the characterized remainder
    #     TODO this should be recursive - the remaining variants could have a one-by-one parametrization
    allele_wo_common = find_corresponding_allele_wo_homozygotes(cursor, variant_ids)
    if allele_wo_common is not None and allele_is_characterized(cursor, allele_wo_common):
        raw = get_params_for_characterized_allele(cursor, allele_wo_common)
        source = AlleleScoringSource(method="characterized_allele_wo_common_variants", allele_id=allele_wo_common,
                                     exp_characterization_ids=exp_characterization_ids(raw))
        return norm_allele_prms_avgd_over_sources(raw, verbose), source

    #     the missing step: take the per-variant minimum of the characterized individual variants
    return parametrize_from_individual_variants(cursor, allele_id, variant_ids, verbose)


def parametrize_genotype(cursor, allele_ids, verbose: bool = False
                         ) -> Tuple[Dict[int, Optional[NormalizedParametrization]], Dict[int, AlleleScoringSource]]:
    """Per-allele normalized parametrization and provenance for a genotype (an iterable of allele ids).
    Keyed by allele id, so a homozygous genotype collapses to a single entry, as before."""
    params_by_allele: Dict[int, Optional[NormalizedParametrization]] = {}
    source_by_allele: Dict[int, AlleleScoringSource] = {}
    for allele_id in allele_ids:
        if allele_id in params_by_allele: continue
        params, source = parametrize_allele(cursor, allele_id, verbose)
        params_by_allele[allele_id] = params
        source_by_allele[allele_id] = source
    return params_by_allele, source_by_allele


def competence_estimate(avg_params: Dict[int, NormalizedParametrization],
                        is_homozygote: bool, assume_dosage_compensation: bool = False,
                        verbose: bool = True) -> float:
    if (is_homozygote and len(avg_params) < 1) or (not is_homozygote and len(avg_params) < 2): return -1

    ability_to_express = {}
    ability_to_transport = {}
    for allele_id, prms in avg_params.items():
        ability_to_express[allele_id] = prms.expression
        # why like this: if the ATPase engine is not working, then the N-Ret-PE_binding site never opens / closes
        # so there is no point in investigating the N-Ret-PE_binding (i.e it should come out as poor even if the
        # binding site is unaffected
        # todo - plot  the ATPase activity vs N-Ret-PE_binding and see how strongly they correlate
        if prms.ATPase_activity < 1:
            ability_to_transport[allele_id] = prms.ATPase_activity \
                                              * prms.ATR_binding_to_ECD2
        else:
            ability_to_transport[allele_id] = prms.N_Ret_PE_binding \
                                              * prms.ATR_binding_to_ECD2

        ability_to_transport[allele_id] *= prms.ATR_binding_to_ECD2

    score = 0
    if assume_dosage_compensation:
        norm = sum(ability_to_transport.values())
        upper_limit = max(ability_to_transport.values())
        if norm > 0:
            fraction_in_membrane = {allele: ability / norm * upper_limit for allele, ability in
                                    ability_to_express.items()}
        else:
            fraction_in_membrane = {allele: 0.0 for allele in ability_to_express.keys()}
    else:
        fraction_in_membrane = {allele: ability * 0.5 for allele, ability in ability_to_express.items()}

    for allele in ability_to_express.keys():
        # there is no evidence for dosage compensation for ABCA4, so we stick 0.5 here
        score += fraction_in_membrane[allele] * ability_to_transport[allele]
        if verbose:
            print(f"\t allele {allele + 1}:  fraction in membrane %.2f  ability to transport: %.2f " %
                  (fraction_in_membrane[allele], ability_to_transport[allele]))
    if is_homozygote: score *= 2  # this mul by 2 is because  the above loop ran only once for a homozygote

    # the multiplication by 2 is so the score runs from 0 to 1 if there is a single functional allele
    # we know that should be the case because Stargardt's disease
    score = min(2 * score, 1)  # anything bigger than 0.5 is enough
    if verbose: print(f"\t score: {score}")
    if verbose: print("==================================================\n")
    return score
