from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from uctreg.common import run_gpt_celltype
from uctreg.configs import L2_REVIEWED_CLUSTER_MAPS, L3_IMMUNE, L3_STROMAL, QC_THRESHOLDS
from uctreg.paths import ProjectPaths, prepare_output, resolve_data_root


def load_stage(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / "scripts" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


merge_stage = load_stage("03_merge_tissue.py", "merge_stage")
final_stage = load_stage("10_finalize_annotations.py", "final_stage")


class PathTests(unittest.TestCase):
    def test_environment_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            previous = os.environ.get("UC_TREG_DATA_ROOT")
            os.environ["UC_TREG_DATA_ROOT"] = directory
            try:
                self.assertEqual(resolve_data_root(None), Path(directory).resolve())
                self.assertEqual(ProjectPaths(Path(directory)).qc_filtered, Path(directory) / "03_qc_filtered")
            finally:
                if previous is None:
                    os.environ.pop("UC_TREG_DATA_ROOT", None)
                else:
                    os.environ["UC_TREG_DATA_ROOT"] = previous

    def test_output_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.txt"
            path.write_text("existing")
            with self.assertRaises(FileExistsError):
                prepare_output(path, overwrite=False)
            self.assertEqual(prepare_output(path, overwrite=True), path)


class PipelineLogicTests(unittest.TestCase):
    def test_eighty_percent_uses_ceiling(self):
        self.assertEqual(merge_stage.minimum_dataset_count(11, 0.8), 9)
        self.assertEqual(merge_stage.minimum_dataset_count(10, 0.8), 8)

    def test_dataset_id(self):
        self.assertEqual(merge_stage.dataset_id("10_GSE296969_qc.h5ad"), "GSE296969")
        self.assertEqual(merge_stage.dataset_id("8_SCP259_qc.h5ad"), "SCP259")

    def test_preserved_qc_exceptions(self):
        self.assertEqual(QC_THRESHOLDS["5_GSE134649"]["max_mito"], 0.2)
        self.assertEqual(QC_THRESHOLDS["10_GSE296969"]["max_genes"], 5000)

    def test_l3_cluster_keys_are_generated(self):
        for config in [*L3_IMMUNE.values(), *L3_STROMAL.values()]:
            suffix = config.get("cluster_suffix", config["stem"])
            generated = {
                f"leiden_{resolution}_{suffix}"
                for resolution in config["resolutions"]
            }
            self.assertIn(config["cluster_key"], generated)

    def test_reviewed_l2_maps_cover_expected_parent_clusters(self):
        self.assertEqual(L2_REVIEWED_CLUSTER_MAPS["other"]["1"], "Glia")
        self.assertEqual(L2_REVIEWED_CLUSTER_MAPS["immune"]["0"], "T_NK")
        for mapping in L2_REVIEWED_CLUSTER_MAPS.values():
            self.assertTrue(mapping)
            self.assertNotIn(None, mapping.values())

    def test_gpt_requires_environment_key_before_import(self):
        previous = os.environ.pop("AGI_API_KEY", None)
        try:
            with self.assertRaises(RuntimeError):
                run_gpt_celltype({}, "colon")
        finally:
            if previous is not None:
                os.environ["AGI_API_KEY"] = previous


@unittest.skipUnless(importlib.util.find_spec("pandas"), "pandas is not installed")
class AnnotationTableTests(unittest.TestCase):
    def test_annotation_coverage(self):
        final_stage.validate_coverage(["a", "b"], ["b", "a"])
        with self.assertRaises(ValueError):
            final_stage.validate_coverage(["a", "b"], ["a", "c"])

    def test_duplicate_annotation_barcodes_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotations.csv"
            path.write_text("cell_id,L3\nc1,Treg\nc1,Treg\n")
            with self.assertRaises(ValueError):
                final_stage.read_annotation_table(path)


@unittest.skipUnless(importlib.util.find_spec("anndata"), "anndata is not installed")
class SyntheticAnnDataTests(unittest.TestCase):
    def test_unique_barcodes(self):
        import anndata as ad
        import numpy as np
        from uctreg.common import validate_unique_obs_names

        adata = ad.AnnData(np.ones((3, 2)))
        adata.obs_names = ["c1", "c2", "c3"]
        validate_unique_obs_names(adata, "synthetic")
        adata.obs_names = ["c1", "c1", "c3"]
        with self.assertRaises(ValueError):
            validate_unique_obs_names(adata, "synthetic")


if __name__ == "__main__":
    unittest.main()
