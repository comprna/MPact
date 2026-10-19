"""Portable CLI regressions using synthetic references and bundled models."""
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from output_columns import bundled_model_path, site_output_columns, variant_output_columns


class SchemaTests(unittest.TestCase):
    def test_compact_and_independent_annotations(self):
        inputs = ['#Chromosome', 'Position', 'Reference', 'Alteration', 'CLNSIG', 'custom']
        default = variant_output_columns(inputs)
        self.assertEqual(default[:len(inputs)], inputs)
        self.assertEqual(len(default), len(set(default)))
        self.assertFalse(any('accessibility' in c or 'conservation' in c for c in default))
        for name in ['delta_zscore', 'delta_p_two_sided', 'non_A_score_note',
                     'mpact_threshold_change', 'mpact_validation_MCC_cutoff']:
            self.assertNotIn(name, default)
        self.assertIn('mpact_delta_alt_minus_ref', default)
        self.assertIn('ref_center_is_A', default)
        self.assertIn('mpact_alt_accessibility_20nt', variant_output_columns(inputs, accessibility=True))
        self.assertNotIn('variant_conservation_score', variant_output_columns(inputs, accessibility=True))
        self.assertIn('variant_conservation_score', variant_output_columns(inputs, conservation=True))
        self.assertNotIn('mpact_alt_accessibility_20nt', variant_output_columns(inputs, conservation=True))
        self.assertNotIn('ref_scan_seq', default)
        self.assertIn('ref_scan_seq', variant_output_columns(inputs, full=True))
        self.assertEqual(len(site_output_columns(['chr', 'start', 'strand'])), 9)

    def test_bundled_model_hashes(self):
        manifest = json.loads((ROOT / 'models/model_manifest.json').read_text())
        self.assertEqual(manifest['default_window_size'], 501)
        for window in [101, 201, 501, 801, 1001]:
            path = Path(bundled_model_path(window))
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                             manifest['models'][str(window)]['sha256'])


@unittest.skipUnless(importlib.util.find_spec('tensorflow') and importlib.util.find_spec('pysam'),
                     'TensorFlow and pysam are required for CLI integration tests')
class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import pysam
        cls.temp = tempfile.TemporaryDirectory(prefix='mpact_compact_')
        cls.folder = Path(cls.temp.name)
        cls.fasta = cls.folder / 'hg38.fa'
        cls.fasta.write_text('>chr1\n' + 'ACGT' * 750 + '\n>chr2\n' + 'ACGT' * 750 + '\n')
        pysam.faidx(str(cls.fasta))
        cls.gtf = cls.folder / 'Homo_sapiens.GRCh38.110.gtf'
        cls.gtf.write_text(
            'chr1\ttest\tgene\t1\t3000\t.\t+\t.\tgene_id "g1"; gene_name "PLUS";\n'
            'chr2\ttest\tgene\t1\t3000\t.\t-\t.\tgene_id "g2"; gene_name "MINUS";\n')
        cls.redi = cls.folder / 'redi.tsv.gz'
        with gzip.open(cls.redi, 'wt') as handle:
            handle.write('Region\tPosition\tStrand\nchr1\t1501\t+\n')
        cls.variants = cls.folder / 'variants.tsv'
        cls.variants.write_text('#Chromosome\tPosition\tReference\tAlteration\tgene_symbols\tCLNSIG\tstrand\tcustom\n'
                                'chr1\t1501\tA\tT\tPLUS\tUncertain_significance\t+\t001\n'
                                'chr1\t1502\tC\tA\tPLUS\tUncertain_significance\t+\tNA\n'
                                'chr2\t1504\tT\tC\tMINUS\tUncertain_significance\t-\t003\n')
        cls.sites = cls.folder / 'sites.tsv'
        cls.sites.write_text('chr\tstart\tstrand\tcustom\nchr1\t1501\t+\tone\n'
                             'chr2\t1504\t-\ttwo\nchr1\t1502\t+\tnonA\nchr1\t99999\t+\tinvalid\n')
        cls.environment = dict(os.environ, CUDA_VISIBLE_DEVICES='-1', TF_CPP_MIN_LOG_LEVEL='3',
                               TF_NUM_INTRAOP_THREADS='2', TF_NUM_INTEROP_THREADS='1',
                               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def call(self, script, args, expected=0):
        answer = subprocess.run([sys.executable, str(ROOT / script)] + list(map(str, args)),
                                cwd=ROOT, env=self.environment, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        self.assertEqual(answer.returncode, expected, answer.stdout)
        return answer.stdout

    def variant_args(self, out):
        return ['--input', self.variants, '--output-tsv', out, '--fasta', self.fasta,
                '--gtf', self.gtf, '--rediportal-gz', self.redi, '--batch-size', 4,
                '--input-chunk-size', 1]

    def site_args(self, out):
        return ['--input', self.sites, '--output-tsv', out, '--fasta', self.fasta,
                '--summary-json', str(out) + '.json', '--summary-tsv', str(out) + '.summary.tsv',
                '--batch-size', 4, '--chunksize', 2]

    def test_variant_default_and_resume(self):
        out = self.folder / 'variant.tsv'
        self.call('score_mpact.py', self.variant_args(out) + ['--keep-raw-temp'])
        with out.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            columns, rows = reader.fieldnames, list(reader)
        self.assertTrue(rows)
        self.assertIn('CLNSIG', columns)
        self.assertIn('strand', columns)  # Original input strand is preserved.
        self.assertNotIn('accessibility_source', columns)
        self.assertNotIn('variant_conservation_score', columns)
        self.assertNotIn('delta_p_two_sided', columns)
        self.assertEqual({r['strand_gencode'] for r in rows}, {'+', '-'})
        self.assertEqual({r['custom'] for r in rows}, {'001', 'NA', '003'})
        for row in rows:
            self.assertAlmostEqual(float(row['mpact_alt_score']) - float(row['mpact_ref_score']),
                                   float(row['mpact_delta_alt_minus_ref']), places=6)
            if row['alt_center_is_A'] == 'False':
                self.assertEqual(float(row['mpact_alt_score']), 0)
        # Preserve raw rows, synthesize a completed checkpoint and repeat finalization.
        checkpoint = Path(str(out) + '.checkpoint.json')
        checkpoint.write_text(json.dumps({'last_input_row': 2, 'candidates_written': len(rows)}))
        self.call('score_mpact.py', self.variant_args(out) + ['--resume', '--keep-raw-temp'])
        with out.open() as handle:
            self.assertEqual(list(csv.DictReader(handle, delimiter='\t')), rows)
        self.call('score_mpact.py', self.variant_args(out) +
                  ['--resume', '--window-size', 101], expected=1)

    def test_site_default_non_A_invalid_and_other_window(self):
        out = self.folder / 'sites.tsv.out'
        self.call('m6A_scorer.py', self.site_args(out))
        with out.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            columns, rows = reader.fieldnames, list(reader)
        self.assertEqual(len(rows), 3)
        self.assertEqual(len(columns), 10)
        self.assertEqual(rows[2]['mpact_score'], '')
        self.assertEqual(rows[2]['scoreable_center_A'], 'False')
        summary = json.loads(Path(str(out) + '.json').read_text())
        self.assertEqual(summary['scoreable_center_A_rows'], 2)
        self.assertEqual(summary['invalid_or_unfetchable_rows'], 1)
        self.assertIsNone(summary['accessibility_source'])
        self.assertIsNone(summary['conservation_source'])
        other = self.folder / 'sites101.tsv'
        self.call('m6A_scorer.py', self.site_args(other) + ['--window-size', 101])
        self.assertEqual(json.loads(Path(str(other) + '.json').read_text())['window_size'], 101)

    @unittest.skipUnless(importlib.util.find_spec('RNA') and importlib.util.find_spec('pyBigWig'),
                         'Optional annotation dependencies not installed')
    def test_annotations_and_parallel_site_folding(self):
        import pyBigWig
        track = self.folder / 'test.bw'
        with pyBigWig.open(str(track), 'w') as handle:
            handle.addHeader([('chr1', 3000), ('chr2', 3000)])
            handle.addEntries(['chr1', 'chr2'], [0, 0], ends=[3000, 3000], values=[0.25, 0.75])
        out = self.folder / 'sites_annotated.tsv'
        self.call('m6A_scorer.py', self.site_args(out) + ['--accessibility',
                  '--accessibility-processes', 2, '--conservation-bigwig', track])
        with out.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(float(rows[0]['site_conservation_score']), .25)
        self.assertTrue(all(0 <= float(r['mpact_ref_accessibility_10nt']) <= 1 for r in rows))
        out = self.folder / 'variant_annotated.tsv'
        self.call('score_mpact.py', self.variant_args(out) + ['--accessibility',
                  '--conservation-bigwig', track])
        with out.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            self.assertIn('a_site_conservation_score', reader.fieldnames)
            self.assertIn('mpact_delta_accessibility_20nt', reader.fieldnames)
            rows = list(reader)
        self.assertTrue(rows)


if __name__ == '__main__':
    unittest.main()
