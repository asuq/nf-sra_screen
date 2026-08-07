# nf-sra_screen 0.4.0 software provenance

Version selection was completed on 6 August 2026. Registry publication and
digest verification were completed on 7 August 2026.

## Scope

This record covers direct pipeline tools, workflow runtimes, custom-image bases,
and explicitly installed scientific libraries. Transitive packages are captured
by immutable image digests and the custom-image inventories. Database/reference
assets, including the existing Sandpiper and GTDB filenames, are unchanged and
deferred to the future `db_prep` work.

- [`software.tsv`](software.tsv) is the direct-version source of truth.
- [`images.tsv`](images.tsv) maps every configured process image to its digest,
  architecture, processes, direct versions, build date, and inventory.
- [`packages/`](packages/) contains `conda list --explicit` and
  `pip freeze --all` output exported from each published custom-image digest.

## Custom builds

All nine custom images were cross-built and loaded locally as `linux/amd64`
images with Docker Buildx. Their common base is
`condaforge/miniforge3:26.3.2-3@sha256:95195a50e20f7929a31a68ea407cefe117423258639efb6b66f0518e4b86f8ee`,
the registry's `linux/amd64` platform digest. New versioned tags were used; no
existing tag was overwritten.

The custom Myloasm image combines Myloasm 0.6.0 with mylotools 2.1.0 so
`final_contig_graph.gfa` can be reconciled with the polished, filtered assembly
before the pipeline publishes `assembly.gfa`.

COMEBin retains PyTorch 1.13.1/CUDA 11.7. VAMB and LorBin GPU retain PyTorch
2.6.0/CUDA 12.4. LorBin uses upstream commit
`e35c65b3a97fe225b04dc0beaddb85fcc4a1af7c`. These compatibility stacks were not
loosened while moving from Mambaforge to Miniforge.

## Validation

Each custom image passed its Dockerfile build assertions and local
version/import/CLI checks before publication. After publication,
[`verify_custom_images.sh`](../../verify_custom_images.sh) pulled all nine by
their tag-and-digest reference and repeated the checks. The mapping image also
passed a tiny Bowtie2/SAMtools/minimap2 read-to-reference smoke test. Published
manifest digests were independently resolved with Docker Buildx and match
`images.tsv`.

[`test_container_provenance.py`](../../../tests/test_container_provenance.py)
rejects mutable or malformed refs, missing digests, undocumented process images,
direct-version mismatches, missing custom inventories, and workflow-runtime
version drift.

GPU validation is limited to image build, imports, package consistency, and
PyTorch CUDA build metadata. No GPU hardware execution or numerical validation
is claimed. A full reference-dependent scientific workflow was not run; database
preparation and corresponding end-to-end validation remain outstanding for
`db_prep`.
