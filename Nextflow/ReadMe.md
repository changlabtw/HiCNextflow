# Troubleshooting Nextflow & Singularity on HPC Environments

**Date**: August 2026
**Author**: Data/Systems Engineer

This document outlines critical troubleshooting steps and best practices derived from deploying a Nextflow-based 3D Genome/Hi-C pipeline on a High-Performance Computing (HPC) cluster managed by Slurm.

## 1. The "Singularity-in-Singularity" Anti-Pattern

### The Problem
Initially, the pipeline was triggered using a containerized Nextflow instance:
```bash
srun singularity exec nextflow.sif nextflow run hic_main.nf ...
```
However, the pipeline logic dictated that individual processes (e.g., `ALIGN_BWA`, `RUN_BIN3C`) needed to run inside their own specific Singularity containers (`bin3c.sif`, `juicer.sif`). This caused a "Singularity-in-Singularity" (nested container) conflict. The nested Nextflow instance lacked access to the host's `singularity` executable and Slurm scheduler commands (`sbatch`), leading to the following error:
`env: 'singularity': No such file or directory` (Exit code 127).

### The Solution
To resolve this, we decoupled the orchestrator from the compute tasks. Nextflow must run **natively on the host node** (as a driver process), while it delegates tasks by launching Singularity containers natively for each process.

## 2. Bootstrapping Java & Nextflow in a Restricted HPC

### The Problem
Running Nextflow directly on the host required a Java Runtime Environment (JDK 11 or higher). However, the HPC lacked a system-wide Java module (`module avail java` returned empty), preventing the standard `curl -s https://get.nextflow.io | bash` installation script from working (throwing a "cannot find java" error).

### The Solution (Rootless Installation)
We bypassed the need for system-level administrator privileges by establishing a user-space environment.

**Method A: Using Conda (Recommended)**
```bash
conda create -n nextflow_env -c conda-forge -c bioconda nextflow openjdk=17
conda activate nextflow_env
```

**Method B: Manual JDK Download**
1. Downloaded a portable OpenJDK tarball (e.g., Adoptium JDK 17) to a local `~/Soft/` directory.
2. Extracted and exported to PATH: `export PATH="/home/user/Soft/jdk-17/bin:$PATH"`
3. Re-ran the Nextflow download script: `curl -s https://get.nextflow.io | bash`

**Fixing Execution Errors:**
When moving the Nextflow binary across machines, executable permissions are often stripped by default Linux security policies. We resolved `permission denied` and `command not found` errors by:
1. Granting execution rights: `chmod +x nextflow`
2. Invoking it via a relative path: `./nextflow -v` (Best practice: Move it to a personal `~/bin` directory and add it to the `$PATH`).

## 3. Nextflow Operation & Best Practices

### Basic Execution
To run a workflow and leverage containerized environments, the execution command looks like this:
```bash
./nextflow run hic_main.nf -profile slurm_cpu -resume
```

### Key Concepts & Optimization
*   **The `-resume` flag**: Always include this flag during development or production. It utilizes Nextflow's caching mechanism, ensuring that if a pipeline fails, it will only restart from the failed step rather than re-running successful, computationally heavy processes (like alignments).
*   **Local vs. Slurm Executors**:
    *   `-profile local`: Nextflow executes all tasks sequentially or concurrently **on the single node** where it was launched. It does not utilize the wider HPC cluster.
    *   `-profile slurm_cpu`: Nextflow acts dynamically as an orchestrator. It translates each task into a separate `sbatch` job and distributes them across multiple HPC nodes based on available resources.
*   **Driver Job Resource Allocation**: When submitting Nextflow to the Slurm queue via an `sbatch` wrapper script (the driver), request **minimal resources** (e.g., 2 CPUs, 4GB RAM). Do not request high core counts (e.g., 32 CPUs) for the wrapper script. Nextflow will independently request the heavy compute resources dynamically through the settings defined in `nextflow.config`.

### Advanced Tip: Resolving Reference Index Errors (Absolute Paths)
When handling reference genomes (e.g., `.fasta`), Nextflow isolates task working directories. Using `path fasta` inside a `process` definition copies/symlinks *only* the `.fasta` file, leaving out essential index files (like BWA's `.bwt`, `.pac`), resulting in `[E::bwa_idx_load_from_disk]` errors. 
**Solution:** Pass reference genome locations as `val fasta` (a string representing the absolute path) so tools read directly from the main database directory where all index files securely reside.
