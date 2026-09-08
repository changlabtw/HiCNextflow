# 3D Genome & Hi-C Analysis Pipeline (Nextflow DSL2)

A high-performance, reproducible, and containerized pipeline for Hi-C data preprocessing and 3D genomic feature analysis built with **Nextflow DSL2**, **Singularity**, and **Docker**.

---

## 📐 Pipeline Workflow Architecture

The pipeline processes raw paired-end FASTQ reads into contact matrices, multi-resolution balanced Cooler files, A/B compartments, Topologically Associating Domains (TADs), and chromatin loops.

```text
FASTQ (Paired-end)
  │
  ▼
[ ALIGN_BWA ] ──> bwa-mem2 alignment (Output: Name-sorted BAM)
  │
  ▼
[ RUN_PAIRTOOLS ] ──> Hi-C pairs parsing, sorting, deduplication & filtering (.valid.pairs.gz)
  │
  ├───► [ BUILD_COOLER ] ──> cooler cload (10 kb) ──► cooler zoomify (.mcool)
  │         │
  │         ▼
  │     [ RUN_COMPARTMENTS ] ──> cooltools A/B compartment analysis (E1 eigenvector)
  │
  └───► [ RUN_JUICER_PRE ] ──> Juicer Tools pre (.hic matrix generation)
            │
            ├───► [ RUN_ARROWHEAD ] (Parallel task) ──► TAD Calling (CPU)
            │
            └───► [ RUN_HICCUPS ]   (Parallel task) ──► Chromatin Loop Calling (★ Requires GPU)
```

---

## 🧩 Pipeline Design & Module Details (`main.nf`)

The workflow is implemented using **Nextflow DSL2** with modularity and resource optimization across key execution steps:

### 1. `ALIGN_BWA`
- **Tool**: `bwa-mem2 mem -5SP`
- **Function**: Optimized paired-end alignment specifically tuned for 5' chimeric Hi-C reads.
- **Optimization**: Dynamically allocates threads between `bwa-mem2` and `samtools sort -n` to directly produce name-sorted BAM files, minimizing disk I/O.
- **Reference Passing**: The reference genome FASTA path is passed as a string (`val fasta`) to ensure `bwa-mem2` has direct access to index files (`.bwt`, `.pac`, `.ann`, etc.) in the database directory without missing symlinks.

### 2. `RUN_PAIRTOOLS`
- **Tools**: `pairtools parse`, `pairtools sort`, `pairtools dedup`, `pairtools select`
- **Function**: Parses alignments into Hi-C contact pairs (`--walks-policy 5unique`), sorts pairs, marks PCR duplicates, and filters for high-quality valid contacts (`UU`, `UR`, `RU`).
- **Output**: Generates `.valid.pairs.gz` and duplicate statistics (`.dedup.stats`).

### 3. `BUILD_COOLER` & `RUN_COMPARTMENTS`
- **Tools**: `cooler`, `cooltools`, `bioframe`
- **Matrix Generation**:
  - `cooler cload`: Builds the base resolution contact matrix (default: `10 kb`).
  - `cooler zoomify`: Generates multi-resolution contact maps (`10kb`, `20kb`, `50kb`, `100kb`, `250kb`, `500kb`, `1Mb`) with Knight-Ruiz (KR) / iterative correction matrix balancing.
- **Compartment Analysis**:
  - Calculates GC fraction across genomic bins using `bioframe`.
  - Performs eigenvector decomposition (`cooltools.eigs_cis`) at `100 kb` resolution to produce A/B compartment profiles (`.compartments.bedgraph`).

### 4. `RUN_JUICER_PRE`, `RUN_ARROWHEAD`, `RUN_HICCUPS`
- **Tool**: `Juicer Tools` (`juicer_tools.jar`)
- **Matrix Generation (`RUN_JUICER_PRE`)**: Converts valid pairs into `.hic` format.
- **Parallel Feature Calling**:
  - **`RUN_ARROWHEAD`**: Identifies Topologically Associating Domains (TADs) using matrix transformation and corner-score metrics (CPU-bound).
  - **`RUN_HICCUPS`**: Detects focal chromatin loops using GPU acceleration.

---

## ⚙️ Execution Profiles & Resource Management (`nextflow.config`)

The pipeline includes preset execution profiles tailored for local workstations and HPC clusters:

| Profile | Target Environment | Container Engine | GPU Acceleration (`HiCCUPS`) |
| :--- | :--- | :--- | :--- |
| `local` | Standalone Linux Workstation | Singularity | `--nv` flag |
| `local_docker` | Local Docker Host | Docker | `--gpus all` |
| `hpc_slurm` | High-Performance Computing (Slurm) | Singularity | `--gres=gpu:1`, `--nv` |

### Detailed Resource Allocations (`hpc_slurm` profile)

| Process | Core Function | CPUs | Memory | Time Limit | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `ALIGN_BWA` | Read Alignment | 32 | 90 GB | 24h | High CPU compute |
| `RUN_PAIRTOOLS` | Pairs Parsing & Dedup | 24 | 60 GB | 12h | Multi-threaded disk I/O |
| `BUILD_COOLER` | Multi-res Matrix (.mcool) | 16 | 60 GB | 8h | Multi-resolution zoomify |
| `RUN_COMPARTMENTS` | A/B Compartment Calling | 8 | 24 GB | 4h | `cooltools` eigenvector analysis |
| `RUN_JUICER_PRE` | `.hic` File Construction | 32 | 120 GB | 12h | Memory-intensive sorting |
| `RUN_ARROWHEAD` | TAD Detection | 8 | 32 GB | 6h | Matrix computation |
| `RUN_HICCUPS` | Chromatin Loop Detection | 8 | 48 GB | 8h | Requires 1x NVIDIA GPU |

---

## 🚀 Quick Start

### Prerequisites
- **Nextflow**: Version 22.10+ (requires Java 11+ or Java 17)
- **Container Engine**: Singularity (Apptainer) or Docker

### Running the Workflow

Execute with Slurm profile:
```bash
nextflow run main.nf \
    -profile hpc_slurm \
    --fastq_dir "/path/to/NGS_DATA/HiC" \
    --reads "/path/to/NGS_DATA/HiC/*_R{1,2}_001.fastq.gz" \
    --fasta "/path/to/reference_FASTA/Homo_sapiens_assembly38.fasta" \
    --chrom_sizes "/path/to/reference_FASTA/chrom.sizes" \
    --outdir "/path/to/output" \
    -resume
```

> **Tip**: Always provide `-resume` to leverage Nextflow's pipeline caching mechanism, allowing interrupted runs to continue without re-executing completed tasks.

---

## 🛠️ Troubleshooting & HPC Best Practices

### 1. Avoid "Singularity-in-Singularity" Conflicts
Do not run the Nextflow driver binary inside a Singularity container when launching processes that themselves spawn Singularity containers. Nextflow must run natively on the host/login node (or via an uncontainerized Slurm driver script) so it can directly invoke `singularity` and `sbatch`.

### 2. HPC Driver Script Resource Sizing
When submitting the Nextflow master script via `sbatch`, allocate minimal resources to the driver job (e.g., 2 CPUs, 4 GB RAM). Nextflow will dynamically schedule individual compute-heavy tasks according to the specifications in `nextflow.config`.

### 3. Reference Index Locality
Ensure reference FASTA paths are passed as string values (`val fasta`) rather than file objects (`path fasta`) inside process definitions when aligning with BWA. This prevents Nextflow from isolating only the `.fasta` file and allows tools to resolve index files located in the parent reference directory.

### 4. Escaping Backslashes in Embedded Python/Bash Code
Nextflow `script:` blocks that use `${}` interpolation are Groovy triple-double-quoted strings (`"""..."""`), which also interpret backslash escape sequences (`\n`, `\t`, etc.) **before** the embedded Python or bash code is ever written out. A literal `\n` intended for Python (e.g. inside an f-string) will instead be consumed by Groovy and converted into a real newline character, producing errors like:
```
SyntaxError: unterminated string literal
```
**Fix**: double every backslash meant for the embedded language — `\n` → `\\n`, `\t` → `\\t`, `\\` → `\\\\` — so Groovy leaves the literal escape sequence intact for Python/bash to interpret at runtime. This is easy to miss because it's silent unless the resulting raw character happens to break syntax (as `\n` does); a stray `\t` or similar can silently produce subtly wrong output without ever raising an error.

### 5. `RUN_PAIRTOOLS` Failing Late With `writer bzf_close: bug encountered`
On long-running samples, `RUN_PAIRTOOLS` can fail near the very end of a multi-hour run with a bare `writer bzf_close: bug encountered` error from htslib's multithreaded bgzf writer, and exit code `140`. This is easy to misdiagnose as OOM, disk-full, or a corrupted input BAM — rule those out first via `sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed`, `df -h` on the work dir, and `samtools quickcheck` on the input BAM.

If those all come back clean but `Elapsed` is long (many hours) relative to a low `MaxRSS`, the likely cause is that `pairtools sort`'s temporary chunk files were written to the NFS-mounted work directory (its default `--tmpdir` is the current working directory, and neither the process nor `nextflow.config` pins it elsewhere) rather than node-local disk. Sustained multi-threaded (`--nproc 24`) read/write I/O over a busy, shared NFS mount for many hours can hit a transient write hiccup — and htslib's threaded bgzf writer surfaces that as this opaque `bug encountered` panic instead of a legible I/O error.

**Fix**: add `scratch true` to `RUN_PAIRTOOLS` so the task's working directory (and therefore `pairtools sort`'s default temp files) is staged on node-local disk instead of NFS, and pin `--tmpdir` explicitly as a second guarantee:
```groovy
withName: 'RUN_PAIRTOOLS' {
    cpus   = 24
    memory = '60 GB'
    time   = '12h'
    scratch = true
    errorStrategy = { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries    = 2
}
```
Also keep `errorStrategy 'retry'` on this process as a safety net — since the failure is environment-driven rather than deterministic, a retry is often enough to clear it even before any config change takes effect.

### 6. `RUN_COMPARTMENTS` Silently Skips the Saddle Plot (`Expected HYPHEN token missing`)
`RUN_COMPARTMENTS` can complete with exit `0` and still be missing `*_saddle.png`/`*_saddle.npz` from its output. Nextflow won't flag this as a failure because those outputs are declared `optional: true`, and the saddle-plot code is wrapped in a bare `try/except` that prints the error to `.command.log` instead of raising. Check that log for a line like:
```
Failed to generate saddle plot: Expected HYPHEN token missing
```

**Root cause**: `Homo_sapiens_assembly38.fasta` (the GATK GRCh38 resource-bundle reference) includes HLA allele contigs whose names contain colons, e.g. `HLA-A*01:01:01:01`. If `chrom.sizes` was derived from the full FASTA, these contig names end up in the cooler file's chromosome list. `cooltools.expected_cis()`/`saddle()` (unlike `eigs_cis()`, which doesn't hit this path) fetch per-chromosome data using UCSC-style region strings, and `cooler`'s region parser splits on `:` expecting `chrom:start-end` — for a name like `HLA-A*01:01:01:01` it finds a coordinate token after the first colon but no following hyphen, and raises `Expected HYPHEN token missing`.

**Fix**: exclude ALT/decoy/HLA contigs from `chrom.sizes` before it reaches `cooler cload`, e.g.:
```bash
grep -E '^chr([0-9]+|X|Y|M)\s' /path/to/reference_FASTA/chrom.sizes > /path/to/reference_FASTA/chrom.sizes.primary
```
and point `--chrom_sizes` at the filtered file. Alternatively, restrict `view_df` to canonical chromosomes inside `RUN_COMPARTMENTS` itself if the full reference is needed elsewhere in the pipeline. Re-run `RUN_COMPARTMENTS` with `-resume` once fixed.

### 7. `RUN_HICCUPS` Fails With `Cannot run program "nvcc": No such file or directory`
`RUN_HICCUPS` reaches the GPU kernel-compilation step and fails with:
```
jcuda.CudaException: Could not prepare PTX for source file '...'
Caused by: java.io.IOException: Cannot run program "nvcc": error=2, No such file or directory
...
GPU/CUDA Installation Not Detected
```
even on a node with a working GPU and driver (a clean `nvidia-smi` showing 0% utilization and free memory doesn't rule this out).

**Root cause**: `juicer_tools hiccups` uses JCuda, which JIT-compiles a `.cu` kernel at runtime by shelling out to `nvcc` (the CUDA *compiler*) rather than shipping a precompiled kernel. Singularity's `--nv` flag only bind-mounts the host's NVIDIA *driver* (`libcuda.so`, `/dev/nvidia*`) into the container so precompiled CUDA binaries can run — it does not supply a compiler. If the container image (`Dockerfile.juicer`) is built from a CUDA `runtime` base image rather than a `devel` one, `nvcc` simply isn't present anywhere in the container, and `--nv` can't fix that. Confirm the host doesn't have it either before assuming a bind-mount will help:
```bash
srun -p <gpu_partition> -w <gpu_node> which nvcc
```

**Fix**: rebuild the Juicer image from a CUDA `devel` base instead of `runtime`, e.g. in `Dockerfile.juicer`:
```dockerfile
# was: FROM nvidia/cuda:11.2.2-runtime-ubuntu20.04
FROM nvidia/cuda:11.2.2-devel-ubuntu20.04
```
Keep the same CUDA version as before unless there's a reason to change it — a `devel` image built against an older CUDA version still runs fine under a newer host driver (CUDA maintains driver forward-compatibility). Add `nvcc --version` to the image's existing build-time sanity checks so a future regression back to a `runtime` base fails at build time instead of hours into a pipeline run. After rebuilding, regenerate the `.sif` (`singularity build ... docker-daemon://...`) and confirm `params.sif3` in your run points at the new file before re-running `RUN_HICCUPS`.

Note this only surfaces once `RUN_HICCUPS`'s output is `optional: true` and `hiccups` is called with `--ignore-sparsity` — otherwise the sparsity check (see the GPU-adjacent `RUN_ARROWHEAD`/`RUN_HICCUPS` sparsity-abort pattern, worth watching for on both processes) exits before ever reaching the GPU kernel step and masks this error entirely.
