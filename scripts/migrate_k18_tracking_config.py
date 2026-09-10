#!/usr/bin/env python3
"""Migrate an explicitly equivalent legacy K18 BC response config to a new file.

This is an offline migration, never a runtime compatibility fallback. Legacy
ideal/scaled/signed-position/reprojected response studies cannot be migrated by
dropping their keys: reproduce them in their original checkout instead.
Equivalence covers response parameters, not legacy ROOT coordinates: regenerate
BcOut hits with chamber-local-v1 geometry before using the migrated config.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

from config_lifecycle import ConfigAuditError, audit_config_bytes, load_policy, parse_config_text


def migrate(text: str, source: str = "legacy config") -> str:
    parsed = parse_config_text(text, source)
    values: dict[str, str] = {}
    for key, value, _ in parsed:
        if key in values:
            raise ConfigAuditError(f"{source}: duplicate key {key}")
        values[key] = value.strip('"')

    def number(key: str) -> float:
        try:
            value = float(values[key])
            if math.isfinite(value):
                return value
        except (KeyError, ValueError):
            pass
        raise ConfigAuditError(f"{source}: missing/nonfinite numeric {key}")

    if number("G4DCSmearResolutionScale") != 1:
        raise ConfigAuditError("only legacy BC scale=1 is equivalent; use the historical checkout for other responses")
    for key in ("G4DCUseTiltedReadout", "G4DCSmearSignedPosition"):
        if key in values and number(key) != 0:
            raise ConfigAuditError(f"{key} changes the response; refusing an automatic physics change")
    if "G4DCUseTiltedReadoutAngleFactor" in values:
        raise ConfigAuditError("historical angle-factor config requires its original checkout")
    for key in ("PK18", "K18BFTAcceptedAbsPdg", "K18BFTPositionSmearSigma"):
        number(key)  # Do not infer missing beam species or an ideal BFT default.
    pdg = number("K18BFTAcceptedAbsPdg")
    if number("PK18") <= 0 or pdg <= 0 or not pdg.is_integer() or pdg > 2147483647:
        raise ConfigAuditError("migration requires positive PK18 and an explicit positive integer beam PDG")
    if number("K18BFTPositionSmearSigma") < 0:
        raise ConfigAuditError("negative BFT sigma was clamped by the old binary; refusing to guess its intended value")
    if "K18TrackingConfigVersion" in values:
        raise ConfigAuditError("expected an unversioned legacy config, not an already versioned one")

    removed = {
        "G4DCSmearResolutionScale", "G4DCUseTiltedReadout",
        "G4DCSmearSignedPosition", "G4DCSmearResolutionScaleSdcIn",
        "G4DCSmearResolutionScaleSdcOut",
    }
    removed_lines = {line for key, _, line in parsed if key in removed}
    result = "K18TrackingConfigVersion: 2\n" + "\n".join(line for i, line in enumerate(text.splitlines(), 1)
                       if i not in removed_lines) + "\n"
    audit_config_bytes(result.encode(), source=source, policy=load_policy(),
                       application="DstK18TrackingGeant4")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True,
                        help="new config path; an existing file is never overwritten")
    args = parser.parse_args()
    try:
        result = migrate(args.config.read_text(), str(args.config))
        with args.output.open("x") as output:
            output.write(result)
    except (OSError, ConfigAuditError) as error:
        parser.exit(1, f"FAIL: {error}\n")
    print(f"Migrated equivalent BC response to {args.output}; "
          "check with bin/DstK18TrackingGeant4 --check-config before running. "
          "Legacy BcOut ROOT must be regenerated with chamber-local-v1 geometry.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
