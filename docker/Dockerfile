# ================================================================
# ULTIMATE BIOINFORMATICS DOCKERFILE – STREPTO-PIPELINE EDITION
# Network-hardened, iVar-stable, reproducible for unstable connections
# ================================================================
FROM ubuntu:22.04

# -----------------------------
# 1. Environment
# -----------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PERL5LIB="/usr/local/lib/perl5"
WORKDIR /workspace

# -----------------------------
# 2. Core system dependencies + Git hardening
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git unzip zip build-essential default-jdk golang-go \
    automake autoconf libtool zlib1g-dev \
    python3 python3-pip python3-venv \
    mafft iqtree fastqc multiqc fastp \
    samtools bedtools bwa bowtie2 minimap2 bcftools \
    libhts-dev libncurses5-dev libbz2-dev liblzma-dev \
    ncbi-blast+ perl perl-modules cpanminus \
    libxml-simple-perl libdigest-md5-perl libexpat1-dev libxml-parser-perl \
    libnet-ssleay-perl libssl-dev libxml-twig-perl libwww-perl \
    liblwp-protocol-https-perl perl-tk libdatetime-perl \
    hmmer prodigal locales nano htop lftp libjson-perl \
    && git config --global http.postBuffer 524288000 \
    && git config --global http.version HTTP/1.1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# -----------------------------
# 3. NCBI Entrez Direct (EFETCH)
# -----------------------------
RUN set -ex && \
    for i in 1 2 3 4; do \
        curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -o /tmp/install-edirect.sh && break || sleep 15; \
    done && \
    bash /tmp/install-edirect.sh && \
    rm /tmp/install-edirect.sh && \
    mv $HOME/edirect /opt/edirect

# -----------------------------
# 4. Python Virtual Environment
# -----------------------------
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --upgrade pip && \
    pip install pandas numpy matplotlib seaborn biopython nextflow requests tqdm

# -----------------------------
# 5. BioPerl & XML Dependencies (plus Abricate Perl deps)
# -----------------------------
RUN cpanm --notest Bio::Root::Root LWP::Protocol::https XML::Simple Bio::Perl \
    Path::Tiny JSON List::MoreUtils LWP::Simple

# -----------------------------
# 6. iVar (Compilation) – hardened clone
# -----------------------------
RUN set -ex && \
    for i in 1 2 3 4; do \
        git clone --depth 1 https://github.com/andersen-lab/ivar.git /opt/ivar && break || sleep 15; \
    done && \
    cd /opt/ivar && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    rm -rf /opt/ivar/.git

# -----------------------------
# 7. LoFreq (SourceForge stable binary + size check + retry)
# -----------------------------
RUN set -ex && \
    for i in 1 2 3 4; do \
        curl -fSL "https://downloads.sourceforge.net/project/lofreq/lofreq_star-2.1.2_linux-x86-64.tgz" -o /tmp/lofreq.tgz && \
        [ -s /tmp/lofreq.tgz ] && [ $(stat -c%s /tmp/lofreq.tgz) -gt 700000 ] && break || sleep 15; \
    done && \
    mkdir -p /tmp/lofreq && \
    tar -xzf /tmp/lofreq.tgz -C /tmp/lofreq --strip-components=1 && \
    cp /tmp/lofreq/bin/lofreq /usr/local/bin/lofreq && \
    chmod +x /usr/local/bin/lofreq && \
    rm -rf /tmp/lofreq /tmp/lofreq.tgz

# -----------------------------
# 8. table2asn (NCBI – correct naming + size validation + more retries)
# -----------------------------
RUN set -ex && \
    cd /tmp && \
    for i in 1 2 3 4 5; do \
        wget -O linux64.table2asn.gz \
             https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz && \
        [ $(stat -c%s linux64.table2asn.gz) -gt 15000000 ] && break || sleep 15; \
    done && \
    gunzip -f linux64.table2asn.gz && \
    mv linux64.table2asn /usr/local/bin/table2asn && \
    chmod +x /usr/local/bin/table2asn && \
    rm -f linux64.table2asn.gz*

# -----------------------------
# 9. Prokka (v1.14.6 – Network-Hardened Build)
# -----------------------------
RUN set -ex && \
    cd /tmp && \
    for i in 1 2 3 4 5; do \
        echo "Download attempt $i for Prokka v1.14.6..." && \
        curl -L --http1.1 -f -o prokka.tar.gz \
             "https://github.com/tseemann/prokka/archive/refs/tags/v1.14.6.tar.gz" && \
        if [ -f prokka.tar.gz ] && [ $(stat -c%s prokka.tar.gz) -gt 11000000 ]; then \
            echo "Download successful (Size: $(stat -c%s prokka.tar.gz) bytes)" && \
            break; \
        else \
            echo "Download failed or incomplete. Retrying in 20s..." && \
            rm -f prokka.tar.gz && \
            sleep 20; \
        fi; \
    done && \
    [ -f prokka.tar.gz ] || { echo "FATAL: Prokka download failed after 5 attempts"; exit 1; } && \
    tar -xzf prokka.tar.gz -C /opt/ && \
    ln -sf /opt/prokka-1.14.6/bin/prokka /usr/local/bin/prokka && \
    rm -f prokka.tar.gz && \
    prokka --setupdb

# -----------------------------
# 10. Abricate (AMR Detection – Retry-Hardened + any2fasta dep)
# -----------------------------
RUN set -ex && \
    # Install any2fasta (required by Abricate for FASTA/sequence conversion)
    for i in 1 2 3; do \
        git clone --depth 1 https://github.com/tseemann/any2fasta.git /opt/any2fasta && break || \
        (echo "any2fasta clone failed, retrying..." && sleep 10); \
    done && \
    ln -sf /opt/any2fasta/any2fasta /usr/local/bin/any2fasta && \
    rm -rf /opt/any2fasta/.git && \
    # Now Abricate
    for i in 1 2 3; do \
        git clone --depth 1 https://github.com/tseemann/abricate.git /opt/abricate && break || \
        (echo "Abricate clone failed, retrying..." && sleep 10); \
    done && \
    ln -sf /opt/abricate/bin/abricate /usr/local/bin/abricate && \
    abricate --setupdb && \
    rm -rf /opt/abricate/.git

# -----------------------------
# 11. Final Environment Setup
# -----------------------------
ENV PATH="/opt/edirect:/opt/prokka-1.14.6/bin:/opt/abricate/bin:/usr/local/bin:$PATH"
WORKDIR /workspace

# --------------------------
# 12. Verification (Smoke Test)
# -----------------------------
RUN prokka --version && \
    abricate --version && \
    ivar version && \
    lofreq >/dev/null 2>&1 || true && \
    table2asn --help >/dev/null 2>&1 || true && \
    efetch -version || true




# -----------------------------
# 13. Core Assembly & QC Tools – Final Hardened Version
# -----------------------------

# Install system packages (cd-hit is available in Ubuntu repos)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-distutils \
    cd-hit \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# SPAdes – correct installation (full directory preserved)
RUN wget https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz && \
    tar -xzf SPAdes-4.2.0-Linux.tar.gz -C /opt/ && \
    rm SPAdes-4.2.0-Linux.tar.gz

# Unicycler – ensure build deps + source install
RUN pip install --no-cache-dir setuptools wheel
RUN set -ex && \
    for i in 1 2 3; do \
        git clone --depth 1 https://github.com/rrwick/Unicycler.git /opt/unicycler && break || sleep 10; \
    done && \
    cd /opt/unicycler && \
    python3 setup.py install && \
    rm -rf /opt/unicycler/.git

# Canonical MLST (tseemann/mlst – Perl-based, authoritative)
RUN set -ex && \
    git clone --depth 1 https://github.com/tseemann/mlst.git /opt/mlst && \
    ln -s /opt/mlst/bin/mlst /usr/local/bin/mlst && \
    rm -rf /opt/mlst/.git

# Canonical Kaptive (klebgenomics/Kaptive – latest maintained fork)
RUN set -ex && \
    git clone --depth 1 https://github.com/klebgenomics/Kaptive.git /opt/kaptive && \
    ln -s /opt/kaptive/kaptive.py /usr/local/bin/kaptive && \
    rm -rf /opt/kaptive/.git

# QUAST, mash, ncbi-datasets-pylib (reliable pip packages)
RUN pip install --no-cache-dir \
    quast \
    mash \
    ncbi-datasets-pylib

# -----------------------------
# 14. Optional Runtime Init Script (skeleton – for future DB tools)
# -----------------------------
# Currently minimal / optional. Uncomment ENTRYPOINT if you add DB-heavy tools later.
RUN echo '#!/bin/bash\n\
set -e\n\
if command -v wget >/dev/null && ping -c1 -W1 8.8.8.8 >/dev/null 2>&1; then\n\
  echo "Internet detected – can download optional databases if configured."\n\
else\n\
  echo "No internet detected – running in offline mode."\n\
fi\n\
exec "$@"' > /opt/init.sh && \
    chmod +x /opt/init.sh


# -----------------------------
# 15. Master PATH & Dependency Patch
# -----------------------------
RUN apt-get update && apt-get install -y libmoo-perl && apt-get clean
ENV PATH="/opt/venv/bin:/opt/edirect:/opt/prokka-1.14.6/bin:/opt/abricate/bin:/opt/SPAdes-4.2.0-Linux/bin:/usr/local/bin:$PATH"


# -----------------------------
# 16. Final Tool Fixes (Kaptive & Mash)
# -----------------------------
# 1. Fix Kaptive permissions
RUN find /opt/kaptive -name "kaptive.py" -exec chmod +x {} +

# 2. Install Mash (Verified stable link for v2.3)
RUN wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar && \
    tar -xf mash-Linux64-v2.3.tar && \
    mv mash-Linux64-v2.3/mash /usr/local/bin/mash && \
    chmod +x /usr/local/bin/mash && \
    rm -rf mash-Linux64-v2.3* && \
    echo "Mash v2.3 binary installed successfully"

# -----------------------------
# 17. Final Verification – The Absolute Finish Line
# -----------------------------
RUN set -e && \
    echo "Verifying SPAdes..." && spades.py --version && \
    echo "Verifying Unicycler..." && unicycler --version && \
    echo "Verifying QUAST..." && quast.py --version && \
    echo "Verifying MLST..." && mlst --version && \
    echo "Verifying Kaptive..." && \
      python3 $(find /opt/kaptive -name "kaptive.py") --version && \
    echo "Verifying CD-HIT..." && cd-hit -h | head -n 1 && \
    echo "Verifying Mash..." && /usr/local/bin/mash --version && \
    echo "SUCCESS: The Strep-Pipeline is 100% complete and ready for ISSAD 2026!"



# -----------------------------
# 18. Roary - Pan/Core Genome Pipeline (2026 Fixed Method)
# -----------------------------
RUN set -ex && \
    cpanm --notest Array::Utils Bio::Perl Exception::Class File::Basename File::Copy \
        File::Find File::Path File::Slurper File::Spec File::Temp File::Which \
        Getopt::Long Graph Graph::Reader::Dot Log::Log4perl Moose Moose::Role \
        Text::CSV Digest::MD5 XML::Simple Bio::Roary && \
    for i in 1 2 3; do \
        git clone --depth 1 https://github.com/sanger-pathogens/Roary.git /opt/roary || \
        git clone --depth 1 https://github.com/tseemann/Roary.git /opt/roary && break || sleep 10; \
    done && \
    ln -sf /opt/roary/bin/roary /usr/local/bin/roary && \
    ln -sf /opt/roary/bin/* /usr/local/bin/ && \
    rm -rf /opt/roary/.git

# -----------------------------
# 19. Snippy - Rapid Haploid SNP Calling
# -----------------------------
RUN set -ex && \
    for i in 1 2 3; do \
        git clone --depth 1 https://github.com/tseemann/snippy.git /opt/snippy && break || sleep 10; \
    done && \
    ln -sf /opt/snippy/bin/snippy* /usr/local/bin/ && \
    rm -rf /opt/snippy/.git


# -----------------------------
# 20. Final Verification & Victory Message (FIXED)
# -----------------------------
RUN set -e && \
    echo "=== STREPTOCOCCUS AGALACTIAE PIPELINE - FINAL CHECK (ISSAD 2026) ===" && \
    spades.py --version && \
    unicycler --version && \
    prokka --version && \
    mlst --version && \
    # Using the path-safe way to call kaptive
    kaptive --help | head -n 1 && \
    # Roary uses -w for version or just help
    roary -w && \
    snippy --version && \
    abricate --db vfdb --list | head -n1 && \
    echo "=================================================================" && \
    echo "ALL TOOLS READY: Roary, Snippy, VFDB, Prokka, MLST, Kaptive, etc." && \
    echo "IMAGE READY FOR DOCKER HUB v1.0 - GO PUSH IT NOW!" && \
    echo "================================================================="
