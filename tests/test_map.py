import json
import tempfile
import unittest
from pathlib import Path

from mrkr.map import map_gene_ids


class MapGeneIdsTest(unittest.TestCase):
    def test_empty_markers_file_is_a_noop(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            markers_path = tmpdir / "markers.json"
            gene_map_path = tmpdir / "gmap.txt"

            markers_path.write_text("[]", encoding="utf-8")
            gene_map_path.write_text("MKI67\tENSG_MKI67\n", encoding="utf-8")

            results = map_gene_ids(markers_path, gene_map_file=gene_map_path, verbose=True)

            self.assertEqual(results, [])

    def test_alias_and_unicode_variants_map_for_human_records(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            markers_path = tmpdir / "markers.json"
            gene_map_path = tmpdir / "gmap.txt"

            markers = [
                {"organism": "homo_sapiens", "feature_name": "KI67", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "PECAM-1", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "PECAM", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "PDGFRΑ", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "IFN-Γ", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "LAG-3", "feature_id": None},
                {"organism": "homo_sapiens", "feature_name": "CD3", "feature_id": None},
                {"organism": "mus_musculus", "feature_name": "KI67", "feature_id": None},
            ]
            markers_path.write_text(json.dumps(markers), encoding="utf-8")
            gene_map_path.write_text(
                "\n".join(
                    [
                        "MKI67\tENSG_MKI67",
                        "PECAM1\tENSG_PECAM1",
                        "PDGFRA\tENSG_PDGFRA",
                        "IFNG\tENSG_IFNG",
                        "LAG3\tENSG_LAG3",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            results = map_gene_ids(markers_path, gene_map_file=gene_map_path)

            self.assertEqual(results[0]["feature_id"], "ENSG_MKI67")
            self.assertEqual(results[1]["feature_id"], "ENSG_PECAM1")
            self.assertEqual(results[2]["feature_id"], "ENSG_PECAM1")
            self.assertEqual(results[3]["feature_id"], "ENSG_PDGFRA")
            self.assertEqual(results[4]["feature_id"], "ENSG_IFNG")
            self.assertEqual(results[5]["feature_id"], "ENSG_LAG3")
            self.assertIsNone(results[6]["feature_id"])
            self.assertIsNone(results[7]["feature_id"])


if __name__ == "__main__":
    unittest.main()
