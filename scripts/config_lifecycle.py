#!/usr/bin/env python3
"""Fail-closed audit for production-style key/value configuration files."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import yaml


REPO = Path(__file__).resolve().parents[1]
DEFAULT_POLICY = REPO / "config/config_lifecycle.yml"
KEY_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?::|=|\s)\s*(\S.*)$")


class ConfigAuditError(RuntimeError):
    pass


@dataclass(frozen=True)
class ConfigAudit:
    path: str
    key_count: int
    sha256: str
    keys: tuple[str, ...]


def load_policy(path: Path = DEFAULT_POLICY) -> dict[str, object]:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise ConfigAuditError(f"cannot read lifecycle policy {path}: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ConfigAuditError(f"unsupported lifecycle policy: {path}")
    return data


def retired_keys(policy: dict[str, object]) -> dict[str, object]:
    value = policy.get("retired_keys", {})
    if not isinstance(value, dict):
        raise ConfigAuditError("retired_keys must be a mapping")
    return value


def parse_config_text(text: str, source: str) -> list[tuple[str, str, int]]:
    lines = text.splitlines()
    parsed: list[tuple[str, str, int]] = []
    for line_number, raw in enumerate(lines, start=1):
        content = raw.split("#", 1)[0].strip()
        if not content:
            continue
        match = KEY_RE.match(content)
        if not match:
            raise ConfigAuditError(
                f"{source}:{line_number}: malformed config line: {raw.rstrip()}"
            )
        parsed.append((match.group(1), match.group(2).strip(), line_number))
    return parsed


def parse_config(path: Path) -> list[tuple[str, str, int]]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ConfigAuditError(f"cannot read config {path}: {exc}") from exc
    return parse_config_text(text, str(path))


def audit_config_bytes(
    content: bytes,
    *,
    source: str,
    policy: dict[str, object],
    max_keys: int | None = None,
    application: str | None = None,
) -> ConfigAudit:
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ConfigAuditError(f"{source}: config is not UTF-8") from exc
    parsed = parse_config_text(text, source)
    seen: dict[str, int] = {}
    retired = retired_keys(policy)
    if application is not None:
        scoped = policy.get("retired_keys_by_application", {})
        if application not in scoped:
            raise ConfigAuditError(f"unknown application policy: {application}")
        retired = {**retired, **scoped[application]}
    for key, _value, line_number in parsed:
        if key in seen:
            raise ConfigAuditError(
                f"{source}:{line_number}: duplicate key {key!r}; "
                f"first defined at line {seen[key]}"
            )
        seen[key] = line_number
        if key in retired:
            detail = retired[key]
            replacement = detail.get("replacement") if isinstance(detail, dict) else detail
            raise ConfigAuditError(
                f"{source}:{line_number}: retired key {key!r}; {replacement}"
            )
    required = policy.get("required_values_by_application", {}).get(application, {})
    values = {key: value.strip('"') for key, value, _ in parsed}
    for key, expected in required.items():
        try:
            valid = float(values[key]) == float(expected)
        except (KeyError, ValueError):
            valid = False
        if not valid:
            raise ConfigAuditError(f"{source}: {application} requires {key}: {expected}; migrate the legacy config explicitly")
    if max_keys is not None and len(seen) > max_keys:
        raise ConfigAuditError(
            f"{source}: {len(seen)} keys exceeds the release budget {max_keys}"
        )
    digest = hashlib.sha256(content).hexdigest()
    return ConfigAudit(source, len(seen), digest, tuple(seen))


def audit_config(
    path: Path,
    *,
    policy: dict[str, object],
    max_keys: int | None = None,
    application: str | None = None,
) -> ConfigAudit:
    try:
        content = path.read_bytes()
    except OSError as exc:
        raise ConfigAuditError(f"cannot read config {path}: {exc}") from exc
    return audit_config_bytes(
        content, source=str(path), policy=policy, max_keys=max_keys,
        application=application,
    )


def expand_paths(values: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    for value in values:
        matches = [Path(item) for item in glob.glob(value)]
        if not matches:
            candidate = Path(value)
            if candidate.exists():
                matches = [candidate]
        paths.extend(matches)
    unique = sorted({path.resolve() for path in paths})
    if not unique:
        raise ConfigAuditError("no configuration files matched")
    not_files = [str(path) for path in unique if not path.is_file()]
    if not_files:
        raise ConfigAuditError("not regular files: " + ", ".join(not_files))
    return unique


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", help="config paths or glob patterns")
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--max-keys", type=int)
    parser.add_argument("--application", help="optional application-specific retired-key policy")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    try:
        policy = load_policy(args.policy)
        results = [
            audit_config(path, policy=policy, max_keys=args.max_keys,
                         application=args.application)
            for path in expand_paths(args.configs)
        ]
    except ConfigAuditError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps([asdict(item) for item in results], indent=2))
    else:
        for item in results:
            print(f"PASS {item.path}: keys={item.key_count} sha256={item.sha256}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
