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
REPOSITORY = r"[a-z0-9.-]+(?::[0-9]+)?/(?:[a-z0-9._-]+/)*[a-z0-9._-]+"
DIGEST_REF = re.compile(rf"{REPOSITORY}@sha256:[0-9a-f]{{64}}")
PROVENANCE_REF = re.compile(
    rf"(?P<repository>{REPOSITORY}):(?P<tag>[\w][\w.-]{{0,127}})"
    r"@(?P<digest>sha256:[0-9a-f]{64})",
    re.ASCII,
)
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
    """Require digest-only references supported by Apptainer's Docker transport."""
    for image in sorted(images):
        assert DIGEST_REF.fullmatch(image), f"mutable or malformed image ref: {image}"


def provenance_runtime_ref(image: str) -> str:
    """Validate a documented versioned reference and return its runtime identity."""
    match = PROVENANCE_REF.fullmatch(image)
    assert match, f"malformed provenance image ref: {image}"
    assert match['tag'].lower() not in FORBIDDEN_TAGS, (
        f"forbidden mutable tag in {image}"
    )
    return f"{match['repository']}@{match['digest']}"


def test_reference_validation() -> None:
    """Protect runtime syntax and preservation of provenance repository/digest."""
    digest = "sha256:" + "a" * 64
    for repository in ("quay.io/example/tool", "localhost:5000/example/tool"):
        runtime = f"{repository}@{digest}"
        tagged = f"{repository}:1.2.3@{digest}"
        assert_immutable({runtime})
        assert provenance_runtime_ref(tagged) == runtime
        invalid_runtime = (
            tagged, repository, f"{repository}:1.2.3",
            runtime[:-1], runtime + "a", runtime.replace("a" * 64, "g" * 64),
            runtime.replace("sha256:", "sha512:"),
        )
        for image in invalid_runtime:
            try:
                assert_immutable({image})
            except AssertionError:
                continue
            raise AssertionError(f"accepted invalid runtime ref: {image}")
        for image in (runtime, tagged[:-1], f"{repository}:latest@{digest}"):
            try:
                provenance_runtime_ref(image)
            except AssertionError:
                continue
            raise AssertionError(f"accepted invalid provenance ref: {image}")


def main() -> int:
    """Validate container configuration against release provenance."""
    test_reference_validation()
    images = configured_images()
    assert_immutable(images)

    image_rows = read_tsv(PROVENANCE_DIR / "images.tsv")
    documented = {provenance_runtime_ref(row["image_ref"]) for row in image_rows}
    assert len(documented) == len(image_rows), "duplicate repository/digest in images.tsv"
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
