# -----------------------------------------------------------------------------
# Dockerfile – Corigami Web & Worker  (Option B: static data on EFS)
# -----------------------------------------------------------------------------
    FROM mambaorg/micromamba:1.5.7-bullseye

    WORKDIR /
    
    # Base env variables
    ENV MPLBACKEND=Agg \
        PYTHONDONTWRITEBYTECODE=1 \
        PIP_NO_CACHE_DIR=1 \
        PATH="/opt/conda/bin:${PATH}" \
        AWS_REGION=us-east-2 \
        USE_AWS=1 
    
    ENV PATH="/usr/local/bin:${PATH}"

    # ── 1. Bio CLI tools (conda env “base”) ───────────────────────────────────────
    COPY environment-bio.yml /tmp/bio.yml
    RUN micromamba install -n base -y -f /tmp/bio.yml \
     && micromamba install -n base -y pip \
     && micromamba clean -afy
    
    # ── 2. System packages ────────────────────────────────────────────────────────
    USER root
    RUN apt-get update && apt-get install -y --no-install-recommends \
            wget \
            ca-certificates \
            build-essential libbz2-dev liblzma-dev \
            libcurl4-gnutls-dev zlib1g-dev \
     && rm -rf /var/lib/apt/lists/*

    # ── 2.5. Install UCSC bigWigToBedGraph ───────────────────────────────────────
    RUN micromamba install -n base -y \
      -c conda-forge \
      -c bioconda \
      ucsc-bigwigtobedgraph \
    && micromamba clean -afy

    # ── 3. Python requirements ────────────────────────────────────────────────────
    COPY requirements.txt .
    RUN pip install --no-cache-dir -r requirements.txt
    
    # ── 4. Application source code ────────────────────────────────────────────────
    COPY . .

    # ensure Python can see your `app/` package
    ENV PYTHONPATH=/:$PYTHONPATH
    
    # ── 5. Runtime command (data already on EFS) ──────────────────────────────────
    EXPOSE 80
    CMD ["gunicorn", "--workers", "2", "--bind", "0.0.0.0:80", "app:create_app()"]
    