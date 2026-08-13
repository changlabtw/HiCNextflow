# Docker images for the Hi-C pipeline

One Dockerfile per tool — matches the "every tool is its own image" design
used by `config/containers.env` and every `scripts/*.sbatch` job.

| Dockerfile | Tool | Produces |
|---|---|---|
| `Dockerfile.hicpro` | HiC-Pro 3.1.0 | `bwa`, `samtools`, `mergeSAM.py`, `mapped_2_hic_fragments.py`, `digest_genome.py`, `hicpro2juicebox.sh` |
| `Dockerfile.hicup` | HiCUP 0.9.2 | `hicup`, `hicup_mapper`, `hicup_digester`, `bowtie2`, `samtools` |
| `Dockerfile.juicer` | Juicer / Juicer Tools 1.9 (CUDA base) | `bwa`, `juicer_tools.jar` (`pre`, `hiccups`, `arrowhead`) |
| `Dockerfile.hicpipe` | Hi-CPIPE (Yaffe-Tanay hicpipe 0.93) | C/Perl bias-correction binaries — **verify build stanza against the tarball's own docs, see comments in the file** |
| `Dockerfile.hickit` | Hickit | `hickit`, `hickit.js` |
| `Dockerfile.bin3c` | bin3C | `bin3C` CLI |
| `Dockerfile.fanc` | FAN-C 0.9.25 | `fanc` CLI + Python module |

## Build everything at once

```bash
cd docker/
./build_all.sh                              # docker build only
./build_all.sh --to-singularity /opt/containers   # also emit .sif files there
```

## Build one image manually

```bash
docker build -t hicpro:3.1.0 -f Dockerfile.hicpro .
singularity build hicpro_3.1.0.sif docker-daemon://hicpro:3.1.0
# or, from a registry instead of the local docker daemon:
docker tag hicpro:3.1.0 your-registry.example.com/hicpro:3.1.0
docker push your-registry.example.com/hicpro:3.1.0
singularity pull hicpro_3.1.0.sif docker://your-registry.example.com/hicpro:3.1.0
```

## Wiring the result into the pipeline

Once you have `.sif` files, point `../config/containers.env` at them:

```bash
export HICPRO_SIF=/opt/containers/hicpro_3.1.0.sif
export HICUP_SIF=/opt/containers/hicup_0.9.2.sif
export JUICER_SIF=/opt/containers/juicer_1.9.sif
export HICPIPE_SIF=/opt/containers/hicpipe_0.93.sif
export HICKIT_SIF=/opt/containers/hickit_latest.sif
export BIN3C_SIF=/opt/containers/bin3c_latest.sif
export FANC_SIF=/opt/containers/fanc_0.9.25.sif
```

## Before trusting these in production

- **Pin real versions.** Several `ARG ..._VERSION` / `ARG ..._COMMIT` values
  default to `master`/`latest` or a specific tag I've set based on current
  public releases — re-check each against the tool's actual releases page
  before building for a production run, since upstream tags and asset URLs
  do shift over time.
- **`Dockerfile.juicer`**: the `juicer_tools.jar` download URL is the most
  likely thing to go stale (aidenlab has moved hosting for this asset more
  than once). If both `wget` sources in the file 404, download the jar
  manually from https://github.com/aidenlab/juicer/wiki/Download and `COPY`
  it in instead — a commented-out `COPY` line is already there.
- **`Dockerfile.hicpipe`**: this is the classic Yaffe & Tanay `hicpipe`
  package (C + Perl, Nature Genetics 2011), distributed as a plain tarball
  with no GitHub releases or bioconda recipe. The build stanza is a
  best-effort reconstruction — open the tarball after downloading and check
  its actual `README`/`INSTALL` before relying on the `make` step as written.
  If your team is actually using a different, more modern "Hi-CPIPE" (e.g.
  a site-internal tool, or the unrelated `ChenFengling/HiCpipe` BL-Hi-C
  wrapper around Juicer+HiC-Pro), swap in that source instead — the two
  tools are not interchangeable despite the similar name.
- **GPU image size**: `Dockerfile.juicer` starts from an NVIDIA CUDA base
  image (~1-2 GB before adding tools) since HiCCUPS needs `--nv` GPU
  passthrough. If you only ever run Arrowhead/`pre` and never HiCCUPS on
  this image, you can switch the base to a plain `eclipse-temurin:8-jdk` or
  similar to shrink it substantially — GPU access is only needed for HiCCUPS.
- Run each image's built-in sanity-check `RUN` lines as your smoke test —
  if `docker build` completes, those checks already passed once, but re-run
  a quick `docker run --rm <image> <tool> --version` after any rebuild.
