from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from config_lifecycle import ConfigAuditError, audit_config_bytes, load_policy
from migrate_k18_tracking_config import migrate

MINIMAL = """K18TrackingConfigVersion: 2
UNPACK: unpack.xml
DIGIT: digit.xml
CMAP: cmap.xml
DCGEO: geometry.dat
K18TM: matrix.dat
USER: user.dat
PK18: 1.4
K18BFTAcceptedAbsPdg: 321
K18BFTPositionSmearSigma: 0.20
"""


class K18TrackingConfigTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (REPO / "tmp").mkdir(exist_ok=True)
        cls.temp = tempfile.TemporaryDirectory(dir=REPO / "tmp")
        cls.root = Path(cls.temp.name)
        source = cls.root / "check.cc"
        source.write_text('''#include "K18Geant4Config.hh"
#include <iostream>
int main(int argc, char** argv) {
  try { k18geant4::ValidateConfigFile(argv[1]); }
  catch(const std::exception& e) { std::cerr << e.what(); return 1; }
  return 0;
}
''')
        cls.probe = cls.root / "check"
        subprocess.run(["g++", "-std=c++17", "-I", str(REPO / "include"),
                        str(source), "-o", str(cls.probe)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def check(self, text, error=None):
        path = self.root / "input.conf"
        path.write_text(text)
        result = subprocess.run([str(self.probe), str(path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0 if error is None else 1, result.stderr)
        if error:
            self.assertIn(error, result.stderr)

    def test_minimal_and_explicit_zero_bft(self):
        self.check(MINIMAL)
        self.check(MINIMAL.replace("0.20", "0"))

    def test_all_scoped_tombstones_rejected_at_runtime(self):
        policy = load_policy()
        retired = {**policy["retired_keys"],
                   **policy["retired_keys_by_application"]["DstK18TrackingGeant4"]}
        for key in retired:
            with self.subTest(key=key):
                self.check(MINIMAL + f"{key}: 0\n", "retired/unsupported")

    def test_duplicate_unknown_and_malformed_rejected(self):
        self.check(MINIMAL + "PK18: 1.2\n", "duplicate")
        self.check(MINIMAL + "K18BFTPositionSmearSigam: 0.2\n", "unknown")
        self.check(MINIMAL + "missing_value\n", "expected")

    def test_species_cannot_fall_back_to_different_bh_and_bft_defaults(self):
        self.check(MINIMAL.replace("K18BFTAcceptedAbsPdg: 321\n", ""), "missing required")
        for value in ("-321", "211.5", "2147483648"):
            self.check(MINIMAL.replace("321", value), "positive integer")

    def test_unversioned_ideal_config_cannot_silently_gain_smearing(self):
        self.check(MINIMAL.replace("K18TrackingConfigVersion: 2\n", ""), "missing required")
        self.check(MINIMAL.replace("ConfigVersion: 2", "ConfigVersion: 1"), "must be 2")

    def test_invalid_response_and_seed_rejected(self):
        for value in ("nan", "inf", "garbage", "0.2mm"):
            self.check(MINIMAL.replace("0.20", value), "finite number")
        self.check(MINIMAL.replace("0.20", "-0.2"), "nonnegative")
        self.check(MINIMAL.replace("K18BFTPositionSmearSigma: 0.20\n", ""), "missing required")
        for value in ("-1", "3.5", "4294967296"):
            self.check(MINIMAL + f"K18HitSmearSeed: {value}\n", "integer in")

    def test_migration_preserves_only_equivalent_response(self):
        legacy = MINIMAL.replace("K18TrackingConfigVersion: 2\n", "") + "G4DCSmearResolutionScale: 1\nG4DCUseTiltedReadout: 0\n"
        self.assertEqual(migrate(legacy), MINIMAL)
        self.check(migrate(legacy))
        for text in (MINIMAL, legacy.replace("Scale: 1", "Scale: 0"),
                     legacy.replace("Readout: 0", "Readout: 1"),
                     legacy.replace("AbsPdg: 321", "AbsPdg: 0"),
                     legacy.replace("Sigma: 0.20", "Sigma: -0.2"),
                     legacy + "G4DCSmearSignedPosition: 1\n"):
            with self.assertRaises(ConfigAuditError):
                migrate(text)

    def test_offline_retirement_is_application_scoped(self):
        content = b"G4DCSmearSignedPosition: 1\n"
        audit_config_bytes(content, source="s2s", policy=load_policy())
        with self.assertRaisesRegex(ConfigAuditError, "retired"):
            audit_config_bytes(content, source="k18", policy=load_policy(),
                               application="DstK18TrackingGeant4")

    def test_dc_smearing_multipliers_are_retired_for_all_applications(self):
        for key in ("G4DCSmearResolutionScale", "G4DCSmearResolutionScaleSdcIn",
                    "G4DCSmearResolutionScaleSdcOut"):
            for application in (None, "DstK18TrackingGeant4"):
                for value in (0, 1, 2):
                    with self.subTest(key=key, application=application, value=value):
                        with self.assertRaisesRegex(ConfigAuditError, "retired"):
                            audit_config_bytes(f"{key}: {value}\n".encode(),
                                               source="smearing", policy=load_policy(),
                                               application=application)

    def test_readout_toggle_is_retired_for_both_systems(self):
        for application in (None, "DstK18TrackingGeant4"):
            for value in (0, 1):
                with self.subTest(application=application, value=value):
                    with self.assertRaisesRegex(ConfigAuditError, "retired"):
                        audit_config_bytes(
                            f"G4DCUseTiltedReadout: {value}\n".encode(),
                            source="readout", policy=load_policy(), application=application)

    def test_public_k18_examples_are_accepted(self):
        configs = sorted((REPO / "runmanager/runlist").glob("*_k18*.conf"))
        self.assertEqual(len(configs), 7)
        configs.append(REPO / "config/k18_tracking_minimal.conf")
        for path in configs:
            with self.subTest(config=path.name):
                text = path.read_text()
                self.check(text)
                audit_config_bytes(text.encode(), source=path.name, policy=load_policy(),
                                   application="DstK18TrackingGeant4")


if __name__ == "__main__":
    unittest.main()
