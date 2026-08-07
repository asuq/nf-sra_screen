#!/usr/bin/env python3
"""Reject mutable, undocumented, or version-inconsistent process images."""

from __future__ import annotations

import csv
import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE = "0.4.0"
PROVENANCE_DIR = ROOT / "docker" / "provenance" / RELEASE
DIGEST_REF = re.compile(r"^.+:[^/@]+@sha256:[0-9a-f]{64}$")
CONTAINER_ASSIGNMENT = re.compile(r'\bcontainer\s*=\s*"([^"]+)"')
FORBIDDEN_TAGS = {"dev", "latest", "main", "master", "snapshot"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated provenance table and fail on missing cells."""
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows, f"empty provenance table: {path}"
    for line_number, row in enumerate(rows, start=2):
        assert None not in row, f"extra columns in {path}:{line_number}"
        missing = [key for key, value in row.items() if not value]
        assert not missing, f"empty cells {missing} in {path}:{line_number}"
    return rows


def configured_images() -> set[str]:
    """Return every literal process image configured by Nextflow."""
    config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    images = set(CONTAINER_ASSIGNMENT.findall(config))
    assert images, "nextflow.config contains no literal process images"
    return images


def assert_immutable(images: set[str]) -> None:
    """Require a versioned tag and a sha256 digest for each image."""
    for image in sorted(images):
        assert DIGEST_REF.fullmatch(image), f"mutable or malformed image ref: {image}"
        tag = image.rsplit("@", maxsplit=1)[0].rsplit(":", maxsplit=1)[1].lower()
        assert tag not in FORBIDDEN_TAGS, f"forbidden mutable tag in {image}"


def main() -> int:
    """Validate container configuration against release provenance."""
    images = configured_images()
    assert_immutable(images)

    image_rows = read_tsv(PROVENANCE_DIR / "images.tsv")
    documented = {row["image_ref"] for row in image_rows}
    assert len(documented) == len(image_rows), "duplicate image_ref in images.tsv"
    assert images == documented, (
        f"configured but undocumented: {sorted(images - documented)}; "
        f"documented but unused: {sorted(documented - images)}"
    )

    software_rows = read_tsv(PROVENANCE_DIR / "software.tsv")
    versions = {row["software"]: row["version"] for row in software_rows}
    assert len(versions) == len(software_rows), "duplicate software in software.tsv"

    used_software: set[str] = set()
    for row in image_rows:
        assert row["architecture"] == "linux/amd64", row["image_ref"]
        for item in row["direct_software"].split(";"):
            name, separator, version = item.partition("=")
            assert separator and name and version, f"malformed direct_software item: {item}"
            assert versions.get(name) == version, (
                f"version mismatch for {name}: image records {version}, "
                f"software.tsv records {versions.get(name)}"
            )
            used_software.add(name)

        if row["kind"] == "custom":
            assert row["build_date"] != "-", f"missing build date: {row['image_ref']}"
            inventory = ROOT / row["inventory"]
            assert inventory.is_file() and inventory.stat().st_size > 0, (
                f"missing package inventory: {inventory}"
            )
        else:
            assert row["inventory"] == "-", row["image_ref"]

    container_targets = {
        row["software"] for row in software_rows if row["scope"] == "container"
    }
    assert container_targets == used_software, (
        f"unrepresented software targets: {sorted(container_targets - used_software)}; "
        f"undocumented software: {sorted(used_software - container_targets)}"
    )

    with (ROOT / "pixi.toml").open("rb") as handle:
        pixi = tomllib.load(handle)
    nextflow = versions["Nextflow"]
    assert pixi["dependencies"]["nextflow"].lstrip("=") == nextflow
    assert f"export NXF_VER={nextflow}" in (ROOT / "run.sh").read_text(encoding="utf-8")
    config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    assert f"version     = '{RELEASE}'" in config
    assert f"nextflowVersion = '>={nextflow}'" in config

    print(
        f"Validated {len(images)} immutable images and "
        f"{len(versions)} direct software targets for {RELEASE}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
