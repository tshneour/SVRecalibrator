FROM condaforge/miniforge3:latest

LABEL description="SVRecalibrator — Structural Variant Breakpoint Refinement Container"

WORKDIR /opt/sv-recalibrator

# Copy environment and source code
COPY environment.yml pyproject.toml collect.py refine.py summarize_homology.py sv_recalibrator.py batch_run.sh README.md ./

# Create conda environment and install package
RUN conda env create -f environment.yml && \
    conda clean -a -y

ENV PATH=/opt/conda/envs/sv-analysis/bin:$PATH

RUN pip install -e .

ENTRYPOINT ["sv-recalibrator"]
CMD ["--help"]
