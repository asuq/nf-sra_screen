#!/usr/bin/env bash

set -euo pipefail

images=(
  "mapping|quay.io/asuq1617/nf-sra_screen_mapping:0.4.0@sha256:c7b08a1057c79e1fa3f3628531ef661c683b4638cd43d6955523b7d4b1b8f0ff"
  "taxa|quay.io/asuq1617/nf-sra_screen_taxa:0.4.0@sha256:6366a34b26b324e1de543b26bb9d91799433ccddea254e0e8e29c32ca964ba44"
  "unicycler|quay.io/asuq1617/nf-sra_screen_unicycler:0.4.0-spades4.3.0@sha256:8b102ba4eb033f9efa7eeff7dc2a0d8aae6027ee6874a8efd20ccd77ce063717"
  "myloasm|quay.io/asuq1617/nf-sra_screen_myloasm:0.4.0-myloasm0.6.0-mylotools2.1.0@sha256:f39e4ca09704cf855db549af22c732dfb8b43b4ed1efe7003b71a6682989bdc2"
  "iseq|quay.io/asuq1617/iseq:1.9.8-sratools3.4.1-r1@sha256:59e96013353dbecddf483aac16591c591b429219b7857a9b277c3bfcb1060ed1"
  "comebin|quay.io/asuq1617/comebin:1.0.4-cuda117-r1@sha256:114a1eb3804e5636b4af07262c078138963b1a2f2ce9a5be73f367dd420b45d6"
  "vamb|quay.io/asuq1617/vamb:5.0.4-torch2.6.0-cuda124-r1@sha256:f102f9aedca1d3b2681d8d1205f1f496fe194067d2d34f2b4fd7a1a107bb2ab2"
  "lorbin-cpu|quay.io/asuq1617/lorbin:0.1.0-cpu-r1@sha256:e07a20bfba278619ae69b21426067922df440efbe240f76860b9333f971f3ff7"
  "lorbin-gpu|quay.io/asuq1617/lorbin:0.1.0-torch2.6.0-cuda124-r1@sha256:743fb052b15366540030df2d64c9d1368dab0b6524a834deb80983084b97a0a6"
)

for entry in "${images[@]}"; do
  IFS='|' read -r name image <<<"${entry}"
  echo "Pulling and verifying ${name}: ${image}"
  docker pull --platform linux/amd64 "${image}" >/dev/null

  case "${name}" in
    mapping)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        test "$(minimap2 --version)" = "2.31-r1302"
        samtools --version | head -1 | grep -Fx "samtools 1.24"
        bowtie2 --version | head -1 | grep -F "version 2.5.5"
        work=$(mktemp -d)
        trap '\''rm -rf "${work}"'\'' EXIT
        cd "${work}"
        printf ">ref\nACGTACGTACGTACGTACGTACGTACGTACGT\n" >ref.fa
        printf "@read\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n" >read.fastq
        bowtie2-build ref.fa index >/dev/null 2>&1
        bowtie2 -x index -U read.fastq -S aligned.sam >/dev/null 2>&1
        samtools view -b aligned.sam >aligned.bam
        samtools quickcheck aligned.bam
        minimap2 ref.fa ref.fa >/dev/null
      '
      ;;
    taxa)
      docker run --platform linux/amd64 --rm "${image}" python -c '
from importlib.metadata import version
expected = {"biopython": "1.88", "numpy": "2.5.1", "openpyxl": "3.1.5", "pandas": "3.0.5", "PyYAML": "6.0.3", "ujson": "5.13.0"}
assert all(version(name) == wanted for name, wanted in expected.items())
import Bio, numpy, openpyxl, pandas, yaml, ujson
'
      ;;
    unicycler)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        unicycler --version | grep -Fx "Unicycler v0.5.1"
        spades.py --version | grep -Fx "SPAdes genome assembler v4.3.0"
        unicycler --help >/dev/null
      '
      ;;
    myloasm)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        set -euo pipefail
        test "$(myloasm --version)" = "myloasm 0.6.0"
        test "$(mylotools --version)" = "2.1.0"
        mylotools annotate-gfa --help >/dev/null
        work=$(mktemp -d)
        trap '\''rm -rf -- "${work}"'\'' EXIT
        cd "${work}"
        printf ">u1ctg_len-4_circular-no_depth-1-1-1_duplicated-no\nACGT\n" >assembly.fasta
        printf "H\tVN:Z:1.0\nS\tu1ctg\t*\tLN:i:3\tDP:f:1.0\nS\tu2ctg\t*\tLN:i:5\tDP:f:2.0\na\tu1ctg,0-3\tread1\nL\tu1ctg\t+\tu2ctg\t+\t0M\n" >final_contig_graph.gfa
        mylotools annotate-gfa \
          --gfa final_contig_graph.gfa \
          --fasta assembly.fasta \
          --output assembly.gfa \
          >/dev/null
        python -c '\''
from pathlib import Path

lines = Path("assembly.gfa").read_text(encoding="utf-8").splitlines()
good = next(line.split("\t") for line in lines if line.startswith("S\tu1ctg\t"))
filtered = next(line.split("\t") for line in lines if line.startswith("S\tu2ctg\t"))
assert good[2] == "ACGT"
assert {"LN:i:4", "LN_GFA:i:3", "FILTERED:Z:GOOD"} <= set(good)
assert filtered[2] == "*"
assert {"LN:i:5", "FILTERED:Z:FAIL"} <= set(filtered)
assert "L\tu1ctg\t+\tu2ctg\t+\t0M" in lines
assert not any(line.startswith("a\t") for line in lines)
'\''
      '
      ;;
    iseq)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        conda list -n iseq --json | python -c '\''import json, sys; packages = {item["name"]: item["version"] for item in json.load(sys.stdin)}; expected = {"iseq": "1.9.8", "sra-tools": "3.4.1", "aria2": "1.37.0", "pigz": "2.8", "pbzip2": "1.1.13"}; assert all(packages.get(name) == version for name, version in expected.items()), packages'\''
        fasterq-dump --version | grep -F "3.4.1"
        iseq -h >/dev/null
      '
      ;;
    comebin)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        conda list -n comebin --json | python -c '\''import json, sys; packages = {item["name"]: item["version"] for item in json.load(sys.stdin)}; assert packages.get("comebin") == "1.0.4", packages'\''
      '
      docker run --platform linux/amd64 --rm "${image}" python -c '
import torch
assert torch.__version__.startswith("1.13.1")
assert torch.version.cuda == "11.7"
'
      docker run --platform linux/amd64 --rm "${image}" bash -lc 'command -v run_comebin.sh >/dev/null'
      ;;
    vamb)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        pip check
        vamb --version | grep -Fx "5.0.4"
        python -c '\''from importlib.metadata import version; import torch, vamb; expected = {"vamb": "5.0.4", "numpy": "1.26.4", "pyhmmer": "0.10.15", "pyrodigal": "3.6.3", "psutil": "5.9.8", "scikit-learn": "1.6.1"}; assert all(version(name) == wanted for name, wanted in expected.items()); assert torch.__version__.startswith("2.6.0") and torch.version.cuda == "12.4"'\''
      '
      ;;
    lorbin-cpu)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        pip check
        python -c '\''from importlib.metadata import version; import lorbin.lorbin, torch; expected = {"LorBin": "0.1.0", "numpy": "1.23.3", "pandas": "2.2.2", "scikit-learn": "1.1.2", "scipy": "1.13.1"}; assert all(version(name) == wanted for name, wanted in expected.items()); assert torch.__version__.startswith("1.11.0") and torch.version.cuda is None'\''
        command -v hmmsearch prodigal bedtools samtools >/dev/null
      '
      ;;
    lorbin-gpu)
      docker run --platform linux/amd64 --rm "${image}" bash -lc '
        pip check
        python -c '\''from importlib.metadata import version; import lorbin.lorbin, torch; expected = {"LorBin": "0.1.0", "numpy": "1.23.3", "pandas": "2.2.2", "scikit-learn": "1.1.2", "scipy": "1.13.1"}; assert all(version(name) == wanted for name, wanted in expected.items()); assert torch.__version__.startswith("2.6.0") and torch.version.cuda == "12.4"'\''
        command -v hmmsearch prodigal bedtools samtools >/dev/null
      '
      ;;
  esac
done

echo "Verified ${#images[@]} custom images by published digest."
