"""
A-to-I editing site filtering utility for MPact.
Loads REDIportal known A-to-I sites and provides overlap checking.
"""

import gzip
from collections import defaultdict


class AtoIFilter:
    """Index known A-to-I editing sites for fast overlap queries."""

    def __init__(self, rediportal_gz_path):
        """
        Load REDIportal hg38 A-to-I sites into memory-indexed structure.
        
        Args:
            rediportal_gz_path: Path to TABLE1_hg38_v3.txt.gz
        """
        self.sites_by_chrom = defaultdict(list)
        self._load(rediportal_gz_path)

    def _load(self, path):
        """Parse REDIportal and index sites by chromosome."""
        with gzip.open(path, "rt", errors="replace") as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            chrom_idx = idx["Region"]
            pos_idx = idx["Position"]
            strand_idx = idx["Strand"]

            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= max(chrom_idx, pos_idx, strand_idx):
                    continue
                chrom = fields[chrom_idx]  # e.g. "chr1"
                try:
                    pos = int(fields[pos_idx])
                except (ValueError, IndexError):
                    continue
                strand = fields[strand_idx]

                # Store 1-based position and strand
                self.sites_by_chrom[chrom].append((pos, strand))

        # Sort each chromosome's sites for binary search
        for chrom in self.sites_by_chrom:
            self.sites_by_chrom[chrom].sort()

    def overlaps_exact(self, chrom, pos, strand=None):
        """
        Check if position exactly overlaps a known A-to-I site.
        
        Args:
            chrom: Chromosome (e.g., "chr1" or "1")
            pos: 1-based genomic position
            strand: Optional strand filter ("+", "-", or None for any)
        
        Returns:
            bool: True if exact overlap found
        """
        if chrom not in self.sites_by_chrom:
            return False
        
        for site_pos, site_strand in self.sites_by_chrom[chrom]:
            if site_pos == pos:
                if strand is None or site_strand == strand:
                    return True
        return False

    def nearest_distance(self, chrom, pos, strand=None):
        """
        Find distance to nearest A-to-I site (0 if exact overlap).
        
        Args:
            chrom: Chromosome
            pos: 1-based position
            strand: Optional strand filter
        
        Returns:
            int: Minimum distance to any A-to-I site, or None if none found
        """
        if chrom not in self.sites_by_chrom:
            return None
        
        min_dist = None
        for site_pos, site_strand in self.sites_by_chrom[chrom]:
            if strand is not None and site_strand != strand:
                continue
            dist = abs(site_pos - pos)
            if min_dist is None or dist < min_dist:
                min_dist = dist
        
        return min_dist

    def nearby_sites_within(self, chrom, pos, radius_nt=5, strand=None):
        """
        Count A-to-I sites within a given radius.
        
        Args:
            chrom: Chromosome
            pos: 1-based position
            radius_nt: Radius in nucleotides
            strand: Optional strand filter
        
        Returns:
            list: List of (site_pos, site_strand, distance) tuples
        """
        if chrom not in self.sites_by_chrom:
            return []
        
        result = []
        for site_pos, site_strand in self.sites_by_chrom[chrom]:
            if strand is not None and site_strand != strand:
                continue
            dist = abs(site_pos - pos)
            if dist <= radius_nt:
                result.append((site_pos, site_strand, dist))
        
        return sorted(result, key=lambda x: x[2])
