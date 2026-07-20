import sqlite3
import tempfile
import unittest
from pathlib import Path

from mrkr.query import parse_query_heuristic, resolve_query


class QueryResolverTest(unittest.TestCase):
    def test_query_uses_aligned_marker_table_gene_labels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            db_path = tmpdir / "llmarkers.sqlite"
            gene_map = tmpdir / "gmap.txt"
            gene_map.write_text("TREM2\tENSG_TREM2\nIL1B\tENSG_IL1B\n", encoding="utf-8")

            with sqlite3.connect(db_path) as conn:
                conn.executescript(
                    """
                    CREATE TABLE papers (
                        paper_id INTEGER PRIMARY KEY,
                        doi TEXT,
                        title TEXT,
                        year INTEGER,
                        license TEXT,
                        abstract TEXT
                    );
                    CREATE TABLE profiles (
                        profile_id INTEGER PRIMARY KEY,
                        paper_id INTEGER NOT NULL,
                        collection TEXT NOT NULL,
                        organism TEXT,
                        group_name TEXT NOT NULL,
                        text_blob TEXT NOT NULL,
                        paper_context_blob TEXT NOT NULL,
                        gene_names_json TEXT NOT NULL,
                        gene_ids_json TEXT NOT NULL,
                        evidence_sentences_json TEXT NOT NULL,
                        n_genes INTEGER NOT NULL,
                        n_gene_ids INTEGER NOT NULL,
                        n_sentences INTEGER NOT NULL
                    );
                    CREATE TABLE markers (
                        marker_id INTEGER PRIMARY KEY,
                        paper_id INTEGER NOT NULL,
                        organism TEXT,
                        group_name TEXT NOT NULL,
                        feature_name TEXT NOT NULL,
                        feature_id TEXT,
                        source_type TEXT NOT NULL,
                        source_rationale TEXT,
                        data_id TEXT
                    );
                    """
                )
                conn.execute(
                    "INSERT INTO papers VALUES (1, '10.0000/example', 'Tumor macrophage study', 2026, NULL, 'tumor context')"
                )
                # Deliberately misalign gene_names_json and gene_ids_json to catch
                # accidental zip-based display logic.
                conn.execute(
                    """
                    INSERT INTO profiles VALUES (
                        10, 1, 'test', 'homo_sapiens', 'MACROPHAGE',
                        'MACROPHAGE TREM2 macrophage marker',
                        'Tumor macrophage study tumor context',
                        '["IL1B", "TREM2"]',
                        '["ENSG_TREM2", "ENSG_IL1B"]',
                        '["TREM2 marks tumor macrophages."]',
                        2, 2, 1
                    )
                    """
                )
                conn.execute(
                    """
                    INSERT INTO markers (
                        paper_id, organism, group_name, feature_name, feature_id,
                        source_type, source_rationale, data_id
                    ) VALUES (
                        1, 'homo_sapiens', 'MACROPHAGE', 'TREM2', 'ENSG_TREM2',
                        'text', 'TREM2 marks tumor macrophages.', NULL
                    )
                    """
                )
                conn.commit()

            query = parse_query_heuristic(
                "TREM2+ macrophages in tumors",
                gene_map_file=gene_map,
            )
            result = resolve_query(db_path, query, top_n=1)

            self.assertEqual(result["query"]["markers"][0]["feature_id"], "ENSG_TREM2")
            self.assertEqual(result["matches"][0]["marker_similarity"]["shared_genes"], ["TREM2"])
            self.assertEqual(result["matches"][0]["marker_similarity"]["query_coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
