"""Required sequence-accessibility and genomic-conservation annotations."""

from __future__ import annotations

import math
import os

import numpy as np


ACCESSIBILITY_LENGTHS = (5, 10, 20)


class AccessibilityCalculator:
    """Calculate ViennaRNA local unpaired probabilities for a centered sequence."""

    def __init__(self, plfold_window=200, maximum_span=150):
        try:
            import RNA
        except ImportError as exc:
            raise ImportError(
                "Accessibility output requires ViennaRNA. Install dependencies with: "
                "pip install -r requirements.txt"
            ) from exc
        self.RNA = RNA
        self.plfold_window = int(plfold_window)
        self.maximum_span = int(maximum_span)

    def score(self, sequence):
        sequence = str(sequence).upper().replace("T", "U")
        if not sequence or len(sequence) % 2 == 0 or set(sequence) - set("ACGUN"):
            raise ValueError("Accessibility requires a non-empty, odd-length A/C/G/T/U/N sequence")
        n = len(sequence)
        window = min(self.plfold_window, n)
        span = min(self.maximum_span, window)
        unpaired = self.RNA.pfl_fold_up(
            sequence, min(max(ACCESSIBILITY_LENGTHS), n), window, span
        )
        center = n // 2
        result = {
            "unpaired_probability_1nt": float(np.clip(unpaired[center + 1][1], 0, 1))
        }
        for length in ACCESSIBILITY_LENGTHS:
            if length > n:
                result[f"accessibility_{length}nt"] = math.nan
                continue
            end = center + 1 + length // 2
            end = min(max(end, length), n)
            result[f"accessibility_{length}nt"] = float(
                np.clip(unpaired[end][length], 0, 1)
            )
        return result


class ConservationTrack:
    """Random-access single-base scores from a phyloP/phastCons bigWig."""

    def __init__(self, path, label=None):
        try:
            import pyBigWig
        except ImportError as exc:
            raise ImportError(
                "Conservation output requires pyBigWig. Install dependencies with: "
                "pip install -r requirements.txt"
            ) from exc
        self.path = os.path.abspath(path)
        self.label = label or os.path.basename(path)
        self.handle = pyBigWig.open(path)
        if self.handle is None:
            raise OSError(f"Could not open conservation bigWig: {path}")
        self.chromosomes = self.handle.chroms()
        self.cache = {}

    def _resolve_chromosome(self, chromosome):
        chromosome = str(chromosome)
        candidates = [chromosome]
        if chromosome.startswith("chr"):
            candidates.append(chromosome[3:])
        else:
            candidates.append("chr" + chromosome)
        return next((value for value in candidates if value in self.chromosomes), None)

    def score(self, chromosome, position1):
        key = (str(chromosome), int(position1))
        if key in self.cache:
            return self.cache[key]
        resolved = self._resolve_chromosome(chromosome)
        if resolved is None or position1 < 1 or position1 > self.chromosomes[resolved]:
            value = math.nan
        else:
            values = self.handle.values(resolved, int(position1) - 1, int(position1))
            value = float(values[0]) if values and values[0] is not None else math.nan
        self.cache[key] = value
        return value

    def close(self):
        self.handle.close()


def annotate_accessibility(frame, calculator, ref_column, alt_column=None):
    """Add centered accessibility fields to a DataFrame in place."""
    suffixes = ("unpaired_probability_1nt", "accessibility_5nt",
                "accessibility_10nt", "accessibility_20nt")
    cache = {}

    def values(sequence):
        sequence = str(sequence)
        if sequence not in cache:
            cache[sequence] = calculator.score(sequence)
        return cache[sequence]

    reference = [values(sequence) for sequence in frame[ref_column]]
    for suffix in suffixes:
        frame[f"mpact_ref_{suffix}"] = [record[suffix] for record in reference]
    if alt_column is not None:
        alternate = [values(sequence) for sequence in frame[alt_column]]
        for suffix in suffixes:
            ref_name = f"mpact_ref_{suffix}"
            alt_name = f"mpact_alt_{suffix}"
            delta_name = f"mpact_delta_{suffix}"
            frame[alt_name] = [record[suffix] for record in alternate]
            frame[delta_name] = frame[alt_name] - frame[ref_name]
    return frame
