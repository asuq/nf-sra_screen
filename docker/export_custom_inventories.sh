#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_dir="${repo_root}/docker/provenance/0.4.0/packages"
mkdir -p "${output_dir}"

images=(
  "nf-sra_screen_mapping-0.4.0|quay.io/asuq1617/nf-sra_screen_mapping:0.4.0@sha256:c7b08a1057c79e1fa3f3628531ef661c683b4638cd43d6955523b7d4b1b8f0ff|mapping"
  "nf-sra_screen_taxa-0.4.0|quay.io/asuq1617/nf-sra_screen_taxa:0.4.0@sha256:6366a34b26b324e1de543b26bb9d91799433ccddea254e0e8e29c32ca964ba44|taxa"
  "nf-sra_screen_unicycler-0.4.0-spades4.3.0|quay.io/asuq1617/nf-sra_screen_unicycler:0.4.0-spades4.3.0@sha256:8b102ba4eb033f9efa7eeff7dc2a0d8aae6027ee6874a8efd20ccd77ce063717|unicycler"
  "iseq-1.9.8-sratools3.4.1-r1|quay.io/asuq1617/iseq:1.9.8-sratools3.4.1-r1@sha256:59e96013353dbecddf483aac16591c591b429219b7857a9b277c3bfcb1060ed1|iseq"
  "comebin-1.0.4-cuda117-r1|quay.io/asuq1617/comebin:1.0.4-cuda117-r1@sha256:114a1eb3804e5636b4af07262c078138963b1a2f2ce9a5be73f367dd420b45d6|comebin"
  "vamb-5.0.4-torch2.6.0-cuda124-r1|quay.io/asuq1617/vamb:5.0.4-torch2.6.0-cuda124-r1@sha256:f102f9aedca1d3b2681d8d1205f1f496fe194067d2d34f2b4fd7a1a107bb2ab2|vamb"
  "lorbin-0.1.0-cpu-r1|quay.io/asuq1617/lorbin:0.1.0-cpu-r1@sha256:e07a20bfba278619ae69b21426067922df440efbe240f76860b9333f971f3ff7|lorbin"
  "lorbin-0.1.0-torch2.6.0-cuda124-r1|quay.io/asuq1617/lorbin:0.1.0-torch2.6.0-cuda124-r1@sha256:743fb052b15366540030df2d64c9d1368dab0b6524a834deb80983084b97a0a6|lorbin"
)

for entry in "${images[@]}"; do
  IFS='|' read -r name image environment <<<"${entry}"
  output="${output_dir}/${name}.txt"
  temporary="${output}.tmp"
  echo "Exporting ${image}"
  docker run --platform linux/amd64 --rm "${image}" \
    bash -lc "printf '%s\n' '# conda explicit inventory' && conda list -n '${environment}' --explicit && printf '%s\n' '# pip freeze --all' && python -m pip freeze --all" \
    >"${temporary}"
  test -s "${temporary}"
  mv "${temporary}" "${output}"
done
