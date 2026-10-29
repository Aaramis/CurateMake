#!/usr/bin/env python3
"""
Improved CoL Taxonomic Matcher
workflow/scripts/CoL_taxonomic_matcher.py

Matches taxonomies against the Catalogue of Life (CoL) API using a modular,
responsibility-separated design.
"""

import sys
import logging
import json
import time
import re
import threading
import requests
from requests.adapters import HTTPAdapter
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Iterator
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from sqlalchemy import text

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)

@dataclass
class TaxonomyRecord:
    """Structured container for taxonomy data."""
    tax_id: int
    kingdom: Optional[str] = None
    phylum: Optional[str] = None
    class_: Optional[str] = None
    order_: Optional[str] = None
    family: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None

@dataclass
class MatchResult:
    """Result of a taxonomy matching attempt."""
    original_tax_id: int
    col_taxonomy: Optional[Dict[str, str]] = None
    col_id: Optional[str] = None  # accepted CoL name-usage id (for reproducibility)
    success: bool = False
    error_message: Optional[str] = None
    # Short reason when the match is defensible but not certain. It becomes
    # seq_tax.status = 'Check' instead of 'Correct', so an ambiguous harmonisation is
    # visible in the database rather than indistinguishable from a clean one.
    review: Optional[str] = None

class CoLAPIClient:
    """ChecklistBank (Catalogue of Life) client using the `nameusage/search` endpoint.

    Rate limiting is THREAD-SAFE (shared lock) and polite: a realistic default rate,
    HTTP 429 handled with Retry-After/backoff, and a descriptive User-Agent. This
    replaces the previous `/match/nameusage` client whose limiter was not thread-safe
    and could not resolve synonyms (see README / SATIVA_modifications style notes).
    """

    HOST = "https://api.checklistbank.org"

    def __init__(self, max_requests_per_second: float = 8.0, timeout: int = 30,
                 contact: str = "https://github.com/Aaramis/CurateMake",
                 dataset: str = "3LR"):
        self.session = requests.Session()
        adapter = HTTPAdapter(pool_connections=20, pool_maxsize=20, max_retries=3, pool_block=True)
        self.session.mount('http://', adapter)
        self.session.mount('https://', adapter)
        self.session.headers.update({'User-Agent': f'CurateMake/1.0 (ITS reference DB curation; +{contact})'})

        # Thread-safe rate limiter. The cap is a safety ceiling against runaway loops;
        # effective throughput is bounded by worker concurrency, not this rate.
        rps = max(0.5, min(float(max_requests_per_second or 150), 500.0))
        self.request_interval = 1.0 / rps
        self._lock = threading.Lock()
        self._last_request = 0.0
        self.timeout = timeout

        # Metrics
        self.api_calls = 0
        self.failed_calls = 0

        # REPRODUCIBILITY: resolve and pin a specific, DOI'd CoL release. `dataset` may be
        # "3LR" (latest COL release, resolved & pinned at run start), a numeric release key,
        # or "3" (the live, moving COL project — NOT reproducible). All matching queries then
        # hit this fixed snapshot; its version + DOI are recorded in provenance.
        self.provenance = self._resolve_dataset(str(dataset))
        self.dataset_key = self.provenance["dataset_key"]

    def _resolve_dataset(self, dataset: str) -> Dict[str, Any]:
        meta = self._get(f"dataset/{dataset}") or {}
        key = meta.get("key")
        if not key:  # resolution failed -> fall back but flag it
            key = dataset if str(dataset).isdigit() else 3
            logger.warning(f"Could not resolve CoL dataset '{dataset}'; falling back to key {key} (provenance unknown)")
        prov = {
            "requested": dataset,
            "dataset_key": key,
            "alias": meta.get("alias"),
            "version": meta.get("version"),
            "issued": meta.get("issued"),
            "doi": meta.get("doi"),
            "title": meta.get("title"),
            "resolved_at": datetime.now().isoformat(),
            "api_host": self.HOST,
        }
        logger.info("CoL source pinned: %s (%s, issued %s) DOI=%s dataset_key=%s",
                    prov.get("alias"), prov.get("version"), prov.get("issued"), prov.get("doi"), key)
        return prov

    def _throttle(self):
        """Serialize + space out requests across all threads."""
        with self._lock:
            wait = self.request_interval - (time.time() - self._last_request)
            if wait > 0:
                time.sleep(wait)
            self._last_request = time.time()

    def _get(self, path: str, params: Optional[Dict[str, str]] = None) -> Optional[Dict]:
        for attempt in range(4):
            self._throttle()
            try:
                r = self.session.get(f"{self.HOST}/{path}", params=params or {}, timeout=self.timeout)
                with self._lock:
                    self.api_calls += 1
                if r.status_code == 200:
                    return r.json()
                if r.status_code == 429:  # rate limited -> respect Retry-After then retry
                    retry = r.headers.get("Retry-After")
                    time.sleep(min(float(retry) if retry else 2 ** attempt, 30) + 0.5)
                    continue
                logger.warning(f"CoL API error {r.status_code} for {path} params: {params}")
                with self._lock:
                    self.failed_calls += 1
                return None
            except requests.RequestException as e:
                logger.warning(f"CoL request failed for {path} {params}: {e}")
                time.sleep(2 ** attempt)
        with self._lock:
            self.failed_calls += 1
        return None

    def search(self, name: str, match_type: str = "EXACT", limit: int = 20) -> List[Dict]:
        """Search the pinned CoL release name index. type=EXACT avoids cross-kingdom epithet
        noise and returns the synonym's accepted usage in the SAME response."""
        d = self._get(f"dataset/{self.dataset_key}/nameusage/search",
                       {"q": name, "type": match_type, "content": "SCIENTIFIC_NAME", "limit": limit})
        return (d or {}).get("result", []) if d else []

    def close(self):
        """Close the HTTP session."""
        self.session.close()

_PROVISIONAL_WORDS = {
    "sp.", "sp", "cf.", "cf", "aff.", "aff", "spp", "spp.", "uncultured", "unidentified",
    "unknown", "unclassified", "environmental", "metagenome", "endophyte", "sample",
    "isolate", "clone", "voucher", "bacterium", "genomosp", "incertae",
}
# CoL rank name -> our TaxonomyRecord field
_RANK_KEY = {"domain": "kingdom", "kingdom": "kingdom", "phylum": "phylum",
             "class": "class_", "order": "order_", "family": "family",
             "genus": "genus", "species": "species"}


def _edit_distance(a: Optional[str], b: Optional[str], cap: int = 3) -> int:
    """Levenshtein distance between two names, short-circuited past `cap`.

    Dependency-free on purpose: this is the only place in the pipeline that needs a
    string metric, and the per-rule conda environments are already large.
    """
    a, b = (a or "").strip().lower(), (b or "").strip().lower()
    if a == b:
        return 0
    if abs(len(a) - len(b)) > cap:
        return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        if min(cur) > cap:
            return cap + 1
        prev = cur
    return prev[-1]


class TaxonomyMatcher:
    """Taxonomy matching against CoL using the EXACT-search strategy.

    `/nameusage/search` with `type=EXACT&content=SCIENTIFIC_NAME` (name only, no context
    in the query) returns the exact usage AND — for a synonym — its accepted name in ONE
    call. Homonyms are disambiguated by filtering the exact results on the expected
    kingdom. Misspellings fall back to a FUZZY search constrained to the same genus.
    Downstream ranks use the same class_/order_ keys the DB layer expects.
    """

    LEVELS = ["species", "genus", "family", "order_", "class_", "phylum", "kingdom"]
    # Maximum edit distance tolerated between a FUZZY query and the name CoL matched.
    MAX_FUZZY_EDITS = 2

    def __init__(self, api_client: CoLAPIClient):
        self.api_client = api_client

    @staticmethod
    def _clean(value: Optional[str]) -> Optional[str]:
        if not value:
            return None
        v = str(value).strip()
        if v.lower() in ('', 'na', 'n/a', 'null', 'none', 'unknown', 'unclassified', 'unidentified'):
            return None
        v = re.sub(r"\s*\(.*", "", v)          # drop author citation, e.g. "(Schaeff.) Pers."
        v = v.replace('_', ' ').strip()
        return v or None

    @staticmethod
    def _is_provisional(name: Optional[str]) -> bool:
        if not name:
            return True
        words = [w.strip(".,;:") for w in name.lower().split()]
        if any(w in _PROVISIONAL_WORDS for w in words):
            return True
        if len(words) >= 2:                     # strain-code epithet, e.g. "sp. F123"
            ep = words[1]
            if ep and sum(c.isdigit() for c in ep) / len(ep) > 0.4:
                return True
        return False

    @staticmethod
    def _kingdom_of(classification: List[Dict]) -> Optional[str]:
        for c in classification or []:
            if c.get("rank") == "kingdom":
                return c.get("name")
        return None

    @staticmethod
    def _has_rank(result: Dict, rank: str, value: str) -> bool:
        return any(c.get("rank") == rank and (c.get("name", "") or "").lower() == value.lower()
                   for c in result.get("classification", []))

    def _narrow(self, cand: List[Dict], rank: str, value: Optional[str]) -> List[Dict]:
        """Soft preference: keep the candidates matching `value` at `rank`, unless none do.

        Soft on purpose. A hard filter on a rank the source may have wrong would turn a
        usable match into no match at all; a preference only breaks ties.
        """
        if not value:
            return cand
        kept = [r for r in cand if self._has_rank(r, rank, value)]
        if not kept:
            return cand

        def live(r):
            return not self._is_any_synonym((r.get("usage", {}) or {}).get("status"))

        # A source lineage can be wrong wholesale, and then this preference makes the
        # source's error win. CoL answers "Leichhardtia" with the accepted Apocynaceae
        # genus, an accepted snail genus, and a synonym of Callitris in Cupressaceae;
        # UNITE files all six of its records under Cupressaceae, so the preference picked
        # the synonym and the accepted genus was gone before anything else could run
        # (reported by Daniel Lundin on UNITE plants). So the preference may choose
        # between LIVING usages, but it must not let a synonym beat an accepted one.
        # Chromis still works: there both candidates are living, an accepted damselfish
        # and a provisionally accepted cichlid, and the source's family still decides.
        if any(live(r) for r in cand) and not any(live(r) for r in kept):
            return cand
        return kept

    def _pick(self, results: List[Dict], expected_kingdom: Optional[str],
              require_genus: Optional[str] = None, expected_rank: Optional[str] = None,
              query: Optional[str] = None, max_edits: Optional[int] = None,
              require_family: Optional[str] = None, prefer_family: Optional[str] = None,
              prefer_order: Optional[str] = None) -> Tuple[Optional[Dict], Optional[str]]:
        """Filter exact/fuzzy results and return the best remaining candidate.

        Guards, applied in order: kingdom (homonyms), genus (typo noise), the rank of the
        result against the rank that was queried, and, for a FUZZY query only, the edit
        distance between what was asked and what CoL matched. An accepted usage wins over
        a synonym among the survivors.
        """
        cand = results
        if expected_kingdom:
            kf = [r for r in cand
                  if (self._kingdom_of(r.get("classification")) or "").lower() == expected_kingdom.lower()]
            if kf:
                cand = kf
        # Within one kingdom the kingdom cannot break a tie, so fall back on the family and
        # then the order the source declares. Reported by Edi Prifti: a Chromis given with
        # family Cichlidae was still attached to Pomacentridae, because the Pomacentridae
        # usage is "accepted" while the cichlid one is only provisionally accepted, and
        # nothing else was consulted.
        cand = self._narrow(cand, "family", prefer_family)
        cand = self._narrow(cand, "order", prefer_order)
        if require_genus:
            cand = [r for r in cand if any(
                c.get("rank") == "genus" and (c.get("name", "").lower() == require_genus.lower())
                for c in r.get("classification", []))]
        if require_family:
            cand = [r for r in cand if any(
                c.get("rank") == "family" and (c.get("name", "").lower() == require_family.lower())
                for c in r.get("classification", []))]
        if expected_rank:
            # Prefer the candidates that sit at the rank actually queried. A rank we do not
            # map (subgenus, subspecies, variety) is kept ONLY when nothing matches exactly,
            # so a subspecies-only match still works while a subgenus homonym stops
            # competing with the genus: CoL answers "Centrophorus" with the shark genus AND
            # the katydid subgenus "Locusta (Centrophorus)", and "Nectamia" with both the
            # genus and the subgenus "Apogon (Nectamia)" (reported by Edi Prifti on fish).
            exact = [r for r in cand if self._result_rank_key(r) == expected_rank]
            cand = exact or [r for r in cand
                             if self._result_rank_key(r) in (None, expected_rank)]
        if max_edits is not None and query:
            cand = [r for r in cand
                    if _edit_distance(query, self._matched_name(r), max_edits) <= max_edits]
        if not cand:
            return None, None
        accepted = [r for r in cand if r.get("usage", {}).get("status") == "accepted"]
        # With nothing accepted, a PLAIN synonym still beats an ambiguous synonym or a
        # misapplied name: it is the only one whose accepted name we are willing to follow,
        # so it is the only one that can give a deeper rank. CoL answers "Anemone patens"
        # with a plain synonym of Pulsatilla halleri and an ambiguous synonym of
        # Pulsatilla patens; taking the ambiguous one lost the species outright.
        plain = [r for r in cand
                 if (r.get("usage", {}) or {}).get("status") == "synonym"]
        chosen = accepted[0] if accepted else (plain[0] if plain else cand[0])

        note = None
        # A synonym is not a competitor: nomenclature has already decided it points
        # elsewhere. What competes with an accepted usage is another LIVING taxon carrying
        # the same name, so the flag decision is taken on the living candidates alone.
        live = [r for r in cand
                if not self._is_any_synonym((r.get("usage", {}) or {}).get("status"))]
        if len(live) == 1:
            settled = (live[0].get("usage", {}) or {}).get("status") == "accepted"
            if not settled:
                # A merely PROVISIONAL usage against an ambiguous synonym of the same name
                # means CoL itself has no settled opinion, so the pick must not look
                # certain. COL26.10 answers "Moolgarda" exactly that way: a provisionally
                # accepted genus in Mugilidae and an ambiguous synonym of Crenimugil, with
                # no accepted usage of the name anywhere. Rule requested by Edi Prifti on
                # the fish test set.
                #
                # Unless the name IS accepted elsewhere in the full response at the same
                # rank. The name is then nomenclaturally valid and what we faced was a
                # homonym, which the family and order guards above resolved on purpose. CoL
                # answers "Chromis" with an accepted damselfish genus (Pomacentridae), a
                # provisionally accepted cichlid genus and an ambiguous synonym of Tilapia:
                # a record declaring Cichlidae wants the cichlid, and needs no flag.
                want_rank = self._result_rank_key(live[0])
                want_name = (self._matched_name(live[0]) or "").lower()
                settled = any(
                    (r.get("usage", {}) or {}).get("status") == "accepted"
                    and (self._matched_name(r) or "").lower() == want_name
                    and self._result_rank_key(r) == want_rank
                    for r in results)
            if settled:
                # One living taxon and synonyms around it is not a guess. CoL answers
                # "Pentacerotidae" with the accepted fish family plus a synonym of the same
                # name pointing at a sea-star family, and "Cyprinidae" with the carp family
                # plus one pointing at a bivalve family. Both are settled.
                return chosen, None
        # Several survivors that disagree on the name they resolve to OR on where they sit:
        # the pick is then a guess, so say so rather than letting it look settled. The
        # lineage half matters for homonyms inside one kingdom, where every candidate
        # carries the SAME name and only the family tells them apart (Chromis is a
        # damselfish and a cichlid).
        def _fam(r):
            for c in r.get("classification", []):
                if c.get("rank") == "family":
                    return c.get("name")
            return None

        signatures = {(self._resolved_name(r), _fam(r)) for r in cand}
        if len(signatures) > 1:
            shown = sorted(f"{n or '?'} ({f or 'famille inconnue'})" for n, f in signatures)
            if not accepted:
                # Nothing is plainly accepted, so the pick really is a judgement call.
                # COL26.10 answers "Moolgarda" with a provisionally accepted genus in
                # Mugilidae and an ambiguous synonym of Crenimugil: either name is
                # defensible, which is exactly why it must not look settled.
                statuses = sorted({(r.get("usage", {}) or {}).get("status") or "?"
                                   for r in cand})
                note = (f"ambiguous: no accepted usage among {len(cand)} candidates "
                        f"({', '.join(statuses)}) -> {shown[:4]}")
            else:
                note = f"ambiguous: {len(cand)} candidates -> {shown[:4]}"
        else:
            usage = chosen.get("usage", {}) or {}
            # Robust to both shapes CoL uses for a synonym it cannot resolve: no `accepted`
            # block at all, or one whose name carries no scientificName.
            if usage.get("status") == "synonym" and not self._resolved_name(chosen):
                # A synonym CoL does not resolve. Falling back on the synonym's own name is
                # the only option, but it is exactly what harmonisation is supposed to fix.
                note = "unresolved synonym: CoL gives no accepted name"
        return chosen, note

    @staticmethod
    def _is_any_synonym(status: Optional[str]) -> bool:
        """True for every CoL status that points at another taxon rather than naming one.

        CoL has "synonym", "ambiguous synonym" and "misapplied name" alongside "accepted"
        and "provisionally accepted". Only the plain synonym is safe to FOLLOW, which is
        why _resolve and _accepted_id still test it exactly. But all three have to be SEEN
        as non-taxa, in two places: an ambiguous synonym must not look like the name it
        shadows when candidates are compared, and it must not count as a rival living taxon
        when deciding whether a pick is a guess. COL26.10 answers "Moolgarda" with a
        provisionally accepted genus and an ambiguous synonym of Crenimugil, both spelled
        Moolgarda in Mugilidae.
        """
        st = (status or "").lower()
        return "synonym" in st or "misapplied" in st

    @classmethod
    def _resolved_name(cls, result: Dict) -> Optional[str]:
        """The name this result points AT, following any flavour of synonym.

        Used only to compare candidates, never to write a taxonomy: an ambiguous synonym
        is precisely one that must not be followed silently.
        """
        usage = result.get("usage", {}) or {}
        name = ((usage.get("accepted") or {}).get("name")
                if cls._is_any_synonym(usage.get("status")) else usage.get("name"))
        return (name or {}).get("scientificName") if isinstance(name, dict) else None

    @staticmethod
    def _binomial_agrees(species: Optional[str], genus: Optional[str]) -> bool:
        """False when the species' own genus token contradicts the genus CoL resolved.

        Stops the source's binomial being written under a different genus: "Rusula azurea"
        must not land under genus "Russula" when only the genus was corrected, and
        "Carangoides ferdau" must not land under genus "Ferdauia" when CoL moved the
        species. A species or genus we cannot compare is accepted.
        """
        if not species or not genus:
            return True
        parts = species.strip().split()
        return bool(parts) and parts[0].lower() == genus.strip().lower()

    @staticmethod
    def _matched_name(result: Dict) -> str:
        """The name CoL actually matched: for a synonym its own name, not the accepted one.

        This is what a FUZZY query is scored against, since the misspelling was written
        against the name as published, not against its accepted replacement.
        """
        name = (result.get("usage", {}) or {}).get("name")
        return (name or {}).get("scientificName", "") if isinstance(name, dict) else ""

    @staticmethod
    def _result_rank_key(result: Dict) -> Optional[str]:
        """The rank key _resolve will file this result's deepest name under.

        This MUST agree with _resolve, because _pick drops candidates whose key is not the
        rank that was queried. An infraspecific accepted name is filed under species there,
        so it has to answer "species" here too: when it answered None, CoL's plain synonym
        of "Anemone patens" (accepted: Pulsatilla halleri subsp. halleri) was dropped from
        a species query and only the ambiguous synonym survived, which is the one whose
        accepted name we refuse to follow. The species was then lost outright.

        None stays the answer for a rank genuinely outside our seven that is not below
        species: such a result is NOT rejected, because its classification still carries
        usable higher ranks.
        """
        usage = result.get("usage", {}) or {}
        if usage.get("status") == "synonym":
            name = (usage.get("accepted") or {}).get("name")
        else:
            name = usage.get("name")
        if not isinstance(name, dict):
            return None
        rank = (name.get("rank") or "").lower()
        if rank in TaxonomyMatcher.INFRASPECIFIC:
            return "species"
        return _RANK_KEY.get(rank)

    @staticmethod
    def _accepted_id(usage: Dict) -> Optional[str]:
        """The accepted name-usage id of a chosen result (for a synonym, the id of its
        accepted usage, not the synonym's own id)."""
        if usage.get("status") == "synonym":
            return (usage.get("accepted") or {}).get("id")
        return usage.get("id")

    # Ranks CoL puts BELOW species. An accepted name at one of these still tells us the
    # species, which is the deepest rank this pipeline stores.
    INFRASPECIFIC = {"subspecies", "variety", "varietas", "form", "forma", "subvariety",
                     "subform", "infraspecific name", "cultivar"}

    def _resolve(self, result: Dict) -> Dict[str, str]:
        """Build a {kingdom..species} dict from a chosen result, resolving synonym->accepted."""
        usage = result.get("usage", {})
        tax: Dict[str, str] = {}
        # higher lineage (for a synonym, CoL's classification is the ACCEPTED placement)
        for c in result.get("classification", []):
            rank = (c.get("rank") or "").lower()
            key = _RANK_KEY.get(rank)
            if key and c.get("name"):
                tax[key] = c["name"]
        # deepest accepted name
        if usage.get("status") == "synonym":
            acc_name = usage.get("accepted", {}).get("name", {})
        else:
            acc_name = usage.get("name", {})
        if isinstance(acc_name, dict):
            sci = acc_name.get("scientificName")
            rank = (acc_name.get("rank") or "").lower()
            key = _RANK_KEY.get(rank)
            if key and sci:
                tax[key] = sci
            elif sci and rank in self.INFRASPECIFIC:
                # CoL often resolves a species synonym to an INFRASPECIFIC accepted name:
                # "Anemone patens" -> "Pulsatilla halleri subsp. halleri". We have no slot
                # below species, and dropping it left the species empty while the
                # classification still supplied the new genus, so the caller then wrote the
                # source's stale binomial back and the row contradicted itself. Keep the
                # species the infraspecific name belongs to.
                parts = sci.split()
                if len(parts) >= 2:
                    tax["species"] = " ".join(parts[:2])
        # Last guard: never return a binomial that contradicts the genus in the same row.
        # It happens when the accepted name is not followed -- an ambiguous synonym, kept
        # on purpose -- while the lineage around it IS the accepted placement:
        # "Senecio aureus" is an ambiguous synonym of both Synotis cappa and Packera aurea,
        # so CoL answers with genus Synotis and the name Senecio aureus. One of the two has
        # to go, and the name we did not resolve is the one we trust least.
        if tax.get("species") and tax.get("genus") and not self._binomial_agrees(
                tax["species"], tax["genus"]):
            tax.pop("species")
        return tax

    def match_taxonomy(self, taxonomy: TaxonomyRecord) -> MatchResult:
        """Match a single taxonomy against CoL with the EXACT-first strategy."""
        try:
            cleaned = {lvl: self._clean(getattr(taxonomy, lvl, None)) for lvl in self.LEVELS}
            exp_kingdom = cleaned.get("kingdom")
            # a provisional species (sp., cf., strain code...) must not be matched as a species
            if cleaned.get("species") and self._is_provisional(cleaned["species"]):
                cleaned["species"] = None

            for lvl in self.LEVELS:                       # deepest -> shallowest
                q = cleaned.get(lvl)
                if not q or self._is_provisional(q):
                    continue
                hit, note = self._pick(self.api_client.search(q, "EXACT"), exp_kingdom,
                                        expected_rank=lvl,
                                        prefer_family=cleaned.get("family"),
                                        prefer_order=cleaned.get("order_"))
                if hit is None and lvl == "species":
                    # Typo in the epithet: the genus is trusted here, so it is the
                    # tightest constraint available.
                    hit, note = self._pick(self.api_client.search(q, "FUZZY"), exp_kingdom,
                                           require_genus=cleaned.get("genus"),
                                           expected_rank=lvl, query=q,
                                           max_edits=self.MAX_FUZZY_EDITS,
                                           prefer_family=cleaned.get("family"),
                                           prefer_order=cleaned.get("order_"))
                elif hit is None and lvl == "genus":
                    # Typo in the GENUS. It cannot be its own constraint: requiring the raw
                    # genus is exactly what made this case unreachable, so a misspelled
                    # genus was never corrected (reported by Edi Prifti). Key the fallback
                    # on the family instead and let the edit distance filter.
                    hit, note = self._pick(self.api_client.search(q, "FUZZY"), exp_kingdom,
                                           require_family=cleaned.get("family"),
                                           expected_rank=lvl, query=q,
                                           max_edits=self.MAX_FUZZY_EDITS,
                                           prefer_order=cleaned.get("order_"))
                    fixed = self._resolve(hit).get("genus") if hit is not None else None
                    if fixed and cleaned.get("species") and fixed.lower() != q.lower():
                        # A corrected genus makes the source's binomial stale. Retry the
                        # species against the corrected genus: that recovers the whole
                        # binomial instead of leaving the epithet behind.
                        deeper, deeper_note = self._pick(
                            self.api_client.search(cleaned["species"], "FUZZY"), exp_kingdom,
                            require_genus=fixed, expected_rank="species",
                            query=cleaned["species"], max_edits=self.MAX_FUZZY_EDITS,
                            prefer_family=cleaned.get("family"))
                        if deeper is not None:
                            hit, lvl, note = deeper, "species", deeper_note
                if hit is not None:
                    tax = self._resolve(hit)
                    # Keep the original ranks MORE SPECIFIC than the one that matched, and
                    # only those. Filling a rank ABOVE the match from the source is what
                    # makes the harmonised layer disagree with itself: for an incertae
                    # sedis family CoL returns no order, and the source's order used to be
                    # reinstated, so the same taxon came out as Perciformes from one source
                    # and Labriformes from another (reported by Edi Prifti on fish eDNA vs
                    # video, 40 of 102 shared species affected). A rank CoL leaves
                    # unresolved must stay empty.
                    for r in self.LEVELS[:self.LEVELS.index(lvl)]:
                        if not cleaned.get(r) or r in tax:
                            continue
                        if r == "species" and not self._binomial_agrees(cleaned[r], tax.get("genus")):
                            continue
                        tax[r] = cleaned[r]
                    if note:
                        logger.warning("tax_id %s: %s (interroge: %r au rang %s)",
                                       taxonomy.tax_id, note, q, lvl)
                    return MatchResult(original_tax_id=taxonomy.tax_id, col_taxonomy=tax,
                                       col_id=self._accepted_id(hit.get("usage", {})),
                                       success=bool(tax), review=note)

            return MatchResult(original_tax_id=taxonomy.tax_id, success=False)

        except Exception as e:
            logger.error(f"Error matching taxonomy {taxonomy.tax_id}: {e}")
            return MatchResult(original_tax_id=taxonomy.tax_id, success=False, error_message=str(e))

def match_taxonomies_parallel(matcher: "TaxonomyMatcher",
                              taxonomies: List[TaxonomyRecord],
                              max_workers: int = 20) -> List[MatchResult]:
    """Match TaxonomyRecords against CoL concurrently.

    THE single implementation of the concurrent CoL query. Used by
    CoLProcessor.process_batch (taxonomy harmonisation, early in the pipeline)
    and by the SATIVA integration prefetch, so both stages talk to
    ChecklistBank the same way. ChecklistBank serves concurrency well (~470
    req/s); a serial stream is what collapses -- measured at one answer per
    25-90s per process once ~90 clade jobs each waited on their own connection.

    Results come back unordered; use MatchResult.original_tax_id to map them
    back to their input.
    """
    if not taxonomies:
        return []

    results: List[MatchResult] = []
    workers = max(1, min(int(max_workers), 32))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_taxonomy = {
            executor.submit(matcher.match_taxonomy, taxonomy): taxonomy
            for taxonomy in taxonomies
        }

        for future in as_completed(future_to_taxonomy):
            taxonomy = future_to_taxonomy[future]
            try:
                results.append(future.result())
            except Exception as e:
                logger.error(f"Future failed for taxonomy {taxonomy.tax_id}: {e}")
                results.append(MatchResult(
                    original_tax_id=taxonomy.tax_id,
                    success=False,
                    error_message=str(e)
                ))

    return results


class CoLProcessor:
    """High-level processor for CoL matching."""

    def __init__(self, db_manager, config: Dict):
        self.db_manager = db_manager
        self.config = config

        # API configuration (polite defaults for the public ChecklistBank API)
        api_config = config.get('CoL_api', {})
        max_requests = api_config.get('max_requests_per_second', 150)
        # Concurrent CoL requests. ChecklistBank handles high concurrency well (~470 req/s,
        # 0 errors), so allow a generous cap; the shared thread-safe limiter is the safety net.
        max_workers = max(1, min(int(api_config.get('max_workers', 20)), 32))
        # CoL release to pin for reproducibility ("3LR" = latest release, resolved at run start)
        dataset = str(api_config.get('dataset', '3LR'))

        self.api_client = CoLAPIClient(max_requests_per_second=max_requests, dataset=dataset)
        self.matcher = TaxonomyMatcher(self.api_client)
        self.max_workers = max_workers

        # Report metadata — records the exact CoL release (version + DOI) used, so the run
        # is reproducible and citable (addresses reviewer note on CoL provenance).
        self.report = {
            "run_summary": {},
            "parameters": {
                "max_requests_per_second": max_requests,
                "max_workers": max_workers,
                "col_source": self.api_client.provenance,
                **config
            },
            "statistics": {}
        }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.api_client.close()

    def _get_table(self, table_name: str) -> str:
        """Get schema-qualified table name for PostgreSQL, unqualified for SQLite."""
        if self.db_manager.db_type == "postgresql":
            schema = self.db_manager.connection_params.get('schema', 'public')
            return f"{schema}.{table_name}"
        return table_name

    def get_pending_taxonomies(self, backbone: str = 'Raw', status: str = 'Pending') -> List[TaxonomyRecord]:
        """Retrieve pending taxonomies that must be processed."""
        taxonomy_table = self._get_table('taxonomy')
        seq_tax_table = self._get_table('seq_tax')

        query = text(f"""
            SELECT DISTINCT t.tax_id, t.kingdom, t.phylum, t.class_, t.order_,
                   t.family, t.genus, t.species
            FROM {taxonomy_table} t
            INNER JOIN {seq_tax_table} st ON t.tax_id = st.tax_id
            WHERE st.backbone = :backbone AND st.status = :status
            ORDER BY t.tax_id
        """)

        def fetch_rows():
            with self.db_manager.get_session() as session:
                result = session.execute(query, {"backbone": backbone, "status": status})
                return list(result)

        rows = self.db_manager._execute_with_retry(
            "col_match_get_pending_taxonomies",
            fetch_rows
        )

        taxonomies = []
        for row in rows:
            taxonomy = TaxonomyRecord(
                tax_id=row[0],
                kingdom=row[1],
                phylum=row[2],
                class_=row[3],
                order_=row[4],
                family=row[5],
                genus=row[6],
                species=row[7]
            )
            taxonomies.append(taxonomy)

        logger.info(f"Retrieved {len(taxonomies)} taxonomies for CoL matching")
        return taxonomies

    def process_batch(self, taxonomies: List[TaxonomyRecord]) -> List[MatchResult]:
        """Process a batch of taxonomies using a thread pool."""
        return match_taxonomies_parallel(self.matcher, taxonomies, self.max_workers)

    def update_database(self, results: List[MatchResult]):
        """Persist successful matches to the database."""
        successful_matches = [r for r in results if r.success and r.col_taxonomy]

        if not successful_matches:
            return 0

        def update_op():
            # Build taxonomy payload once to avoid per-row get_or_create calls
            taxonomy_data = []
            taxonomy_tuple_by_original = {}
            # Ambiguity travels by original tax_id: the loop that writes seq_tax iterates a
            # SQLAlchemy result set, not these MatchResults, so it cannot read .review off
            # the row it is inserting.
            review_by_original = {r.original_tax_id: r.review for r in successful_matches}
            for result in successful_matches:
                col_taxonomy = result.col_taxonomy or {}
                tax_values = {
                    'kingdom': col_taxonomy.get('kingdom'),
                    'phylum': col_taxonomy.get('phylum'),
                    'class_': col_taxonomy.get('class_'),
                    'order_': col_taxonomy.get('order_'),
                    'family': col_taxonomy.get('family'),
                    'genus': col_taxonomy.get('genus'),
                    'species': col_taxonomy.get('species')
                }
                taxonomy_data.append(tax_values)
                taxonomy_tuple_by_original[result.original_tax_id] = (
                    tax_values['kingdom'],
                    tax_values['phylum'],
                    tax_values['class_'],
                    tax_values['order_'],
                    tax_values['family'],
                    tax_values['genus'],
                    tax_values['species']
                )

            tax_id_map_by_tuple = self.db_manager.bulk_insert_taxonomies(taxonomy_data)
            tax_id_map = {
                original_tax_id: tax_id_map_by_tuple.get(tax_tuple)
                for original_tax_id, tax_tuple in taxonomy_tuple_by_original.items()
            }

            # Update existing relations
            original_ids = [tid for tid in tax_id_map.keys() if tid]
            if not original_ids:
                return 0

            seq_tax_table = self._get_table('seq_tax')

            with self.db_manager.get_session() as session:
                if self.db_manager.db_type == "postgresql":
                    result = session.execute(text(f"""
                        SELECT seq_id, tax_id FROM {seq_tax_table}
                        WHERE tax_id = ANY(:tax_ids) AND backbone = 'Raw'
                    """), {"tax_ids": original_ids})
                else:
                    placeholders = ','.join([f':id{i}' for i in range(len(original_ids))])
                    params = {f'id{i}': tax_id for i, tax_id in enumerate(original_ids)}
                    result = session.execute(text(f"""
                        SELECT seq_id, tax_id FROM {seq_tax_table}
                        WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
                    """), params)

                # Prepare CoL backbone insertions
                col_inserts = []
                for seq_id, original_tax_id in result:
                    new_tax_id = tax_id_map.get(original_tax_id)
                    if new_tax_id:
                        col_inserts.append({
                            'seq_id': seq_id,
                            'tax_id': new_tax_id,
                            'backbone': 'CoL',
                            # 'Check' when the harmonisation is defensible but not certain,
                            # so an ambiguous row is distinguishable from a clean one.
                            'status': ('Check' if review_by_original.get(original_tax_id)
                                       else 'Correct')
                        })

                # Bulk insert seq_tax relations
                if col_inserts:
                    self.db_manager.bulk_insert_seq_tax(col_inserts)

                # Record CoL provenance (accepted usage id + pinned release) on each matched
                # taxonomy row. This also stamps rows that already existed from the Raw
                # ingestion (shared lineage): a plain UPDATE, so it covers both cases.
                taxonomy_table = self._get_table('taxonomy')
                dataset_key = self.api_client.dataset_key
                seen_tid = set()
                col_prov_updates = []
                for result in successful_matches:
                    new_tid = tax_id_map.get(result.original_tax_id)
                    if new_tid and result.col_id and new_tid not in seen_tid:
                        seen_tid.add(new_tid)
                        col_prov_updates.append(
                            {"tid": new_tid, "cid": result.col_id, "dk": dataset_key})
                if col_prov_updates:
                    session.execute(text(f"""
                        UPDATE {taxonomy_table}
                        SET col_id = :cid, col_dataset_key = :dk
                        WHERE tax_id = :tid
                    """), col_prov_updates)

                # Mark legacy relations as correct in one update
                if self.db_manager.db_type == "postgresql":
                    session.execute(text(f"""
                        UPDATE {seq_tax_table} SET status = 'Correct'
                        WHERE tax_id = ANY(:tax_ids) AND backbone = 'Raw'
                    """), {"tax_ids": original_ids})
                else:
                    placeholders = ','.join([f':id{i}' for i in range(len(original_ids))])
                    params = {f'id{i}': tax_id for i, tax_id in enumerate(original_ids)}
                    session.execute(text(f"""
                        UPDATE {seq_tax_table} SET status = 'Correct'
                        WHERE tax_id IN ({placeholders}) AND backbone = 'Raw'
                    """), params)

                session.commit()

            return len(successful_matches)

        try:
            return self.db_manager._execute_with_retry(
                "col_match_update_database",
                update_op
            )
        except Exception as e:
            logger.error(f"Database update failed: {e}")
            raise

    def mark_failed_taxonomies(self):
        """Mark taxonomies that never matched as Failed."""
        seq_tax_table = self._get_table('seq_tax')
        def update_op():
            with self.db_manager.get_session() as session:
                result = session.execute(text(f"""
                    UPDATE {seq_tax_table}
                    SET status = 'Failed'
                    WHERE backbone = 'Raw' AND status = 'Pending'
                """))
                failed_count = result.rowcount
                session.commit()
                return failed_count

        failed_count = self.db_manager._execute_with_retry(
            "col_match_mark_failed",
            update_op
        )
        logger.info(f"Marked {failed_count} taxonomies as Failed")
        return failed_count

    def remove_duplicate_taxonomies(self):
        """Remove duplicate taxonomies from the taxonomy table."""
        taxonomy_table = self._get_table('taxonomy')
        seq_tax_table = self._get_table('seq_tax')

        # Adjust expressions per database type
        if self.db_manager.db_type == "postgresql":
            group_concat_expr = "STRING_AGG(CAST(tax_id AS TEXT), ',')"
        else:
            group_concat_expr = "GROUP_CONCAT(tax_id)"

        query = text(f"""
            SELECT kingdom, phylum, class_, order_, family, genus, species,
                    MIN(tax_id) as keep_id,
                    {group_concat_expr} as all_ids
            FROM {taxonomy_table}
            WHERE parent_tax_id IS NOT NULL
            GROUP BY kingdom, phylum, class_, order_, family, genus, species
            HAVING COUNT(*) > 1
        """)

        def dedupe_op():
            with self.db_manager.get_session() as session:
                result = session.execute(query)
                duplicates = result.fetchall()

                if not duplicates:
                    return 0

                removed_count = 0
                for dup in duplicates:
                    keep_id = dup[7]  # keep_id (UUID string)
                    all_ids = [x.strip() for x in dup[8].split(',')]  # all_ids (UUID strings)
                    remove_ids = [x for x in all_ids if x != keep_id]

                    if remove_ids:
                        # Update seq_tax entries
                        if self.db_manager.db_type == "postgresql":
                            session.execute(text(f"""
                                UPDATE {seq_tax_table} SET tax_id = :keep_id
                                WHERE tax_id = ANY(:remove_ids)
                            """), {"keep_id": keep_id, "remove_ids": remove_ids})

                            session.execute(text(f"""
                                DELETE FROM {taxonomy_table} WHERE tax_id = ANY(:remove_ids)
                            """), {"remove_ids": remove_ids})
                        else:
                            placeholders = ','.join([f':id{i}' for i in range(len(remove_ids))])
                            params = {f'id{i}': rid for i, rid in enumerate(remove_ids)}
                            params['keep_id'] = keep_id

                            session.execute(text(f"""
                                UPDATE {seq_tax_table} SET tax_id = :keep_id
                                WHERE tax_id IN ({placeholders})
                            """), params)

                            session.execute(text(f"""
                                DELETE FROM {taxonomy_table} WHERE tax_id IN ({placeholders})
                            """), {f'id{i}': rid for i, rid in enumerate(remove_ids)})

                        removed_count += len(remove_ids)

                session.commit()
                return removed_count

        removed_count = self.db_manager._execute_with_retry(
            "col_match_remove_duplicates",
            dedupe_op
        )
        if removed_count == 0:
            logger.info("No duplicate taxonomies found")
        else:
            logger.info(f"Removed {removed_count} duplicate taxonomies")
        return removed_count

    def run(self, batch_size: int = 100):
        """Execute the main processing loop."""
        start_time = time.time()
        self.report["run_summary"]["start_time"] = datetime.now().isoformat()
        
        logger.info("Starting CoL taxonomic matching...")
        
        try:
            # Retrieve taxonomies
            taxonomies = self.get_pending_taxonomies()
            total_taxonomies = len(taxonomies)
            
            if not taxonomies:
                logger.info("No taxonomies found for processing")
                self._finalize_report(start_time, 0, 0, 0)
                return
            
            # Batch processing
            total_matched = 0
            for i in range(0, len(taxonomies), batch_size):
                batch = taxonomies[i:i + batch_size]
                results = self.process_batch(batch)
                matched_in_batch = self.update_database(results)
                total_matched += matched_in_batch
                
                progress = (i + len(batch)) / len(taxonomies) * 100
                logger.info(f"Progress: {progress:.1f}% - Batch matched: {matched_in_batch}/{len(batch)}")
            
            # Final cleanup
            duplicates_removed = self.remove_duplicate_taxonomies()
            failed_marked = self.mark_failed_taxonomies()
            
            self._finalize_report(start_time, total_taxonomies, total_matched, failed_marked)
            
            logger.info(f"CoL matching completed successfully")
            logger.info(f"Matched: {total_matched}, Failed: {failed_marked}")
            
        except Exception as e:
            logger.error(f"Error in CoL processing: {e}", exc_info=True)
            self.report["run_summary"]["status"] = "Error"
            self.report["run_summary"]["error_message"] = str(e)
            raise

    def _finalize_report(self, start_time: float, total: int, matched: int, failed: int):
        """Finalize and persist statistics in the in-memory report."""
        end_time = time.time()
        processing_time = end_time - start_time
        
        self.report["run_summary"].update({
            "end_time": datetime.now().isoformat(),
            "total_duration_seconds": round(processing_time, 2),
            "status": "Success"
        })
        
        self.report["statistics"].update({
            "total_taxonomies": total,
            "matched_taxonomies": matched,
            "failed_taxonomies": failed,
            "api_calls": self.api_client.api_calls,
            "failed_api_calls": self.api_client.failed_calls,
            "success_rate": (matched / total * 100) if total > 0 else 0,
            "avg_time_per_call": processing_time / max(self.api_client.api_calls, 1)
        })

    def write_report(self, output_path: str):
        """Write the JSON report to disk."""
        logger.info(f"Writing report to: {output_path}")
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with output_file.open('w', encoding='utf-8') as f:
            json.dump(self.report, f, indent=2, ensure_ascii=False)

def main():
    """Script entry point."""
    try:
        # Snakemake or standalone configuration
        if 'snakemake' in globals():
            db_path = snakemake.params.db
            batch_size = snakemake.params.batch_size
            config = snakemake.config
            output_report = snakemake.output[0]
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log

            # Logging configuration for Snakemake
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            logger.warning("Running in standalone test mode")
            db_path = "test_output.sqlite"
            config = {
                'database': {'path': db_path},
                'CoL_api': {
                    'max_requests_per_second': 10,
                    'max_workers': 4
                }
            }
            batch_size=100
            output_report = "CoL_matching_report.json"

        # Import and initialize the database manager
        sys.path.append(str(Path(__file__).parent / "utils"))
        from db_manager import get_database_manager

        db_manager = get_database_manager(config)

        # Main processing
        with CoLProcessor(db_manager, config) as processor:
            processor.run(batch_size=batch_size)
            processor.write_report(output_report)

            # Display aggregated statistics
            stats = processor.report["statistics"]
            logger.info("\n" + "="*50)
            logger.info("CoL MATCHING SUMMARY")
            logger.info("="*50)
            logger.info(f"Total taxonomies processed: {stats['total_taxonomies']}")
            logger.info(f"Successfully matched: {stats['matched_taxonomies']}")
            logger.info(f"Failed to match: {stats['failed_taxonomies']}")
            logger.info(f"Success rate: {stats['success_rate']:.1f}%")
            logger.info(f"API calls made: {stats['api_calls']}")
            logger.info("="*50)

    except Exception as e:
        logger.critical(f"Fatal error in CoL matching: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
