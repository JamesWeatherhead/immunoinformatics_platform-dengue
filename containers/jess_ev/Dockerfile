FROM continuumio/anaconda3

COPY src /

# Gurobi runtime configuration. Credentials must be supplied at runtime
# (NOT baked into the image). Set via -e WLSACCESSID=... -e WLSSECRET=...
# -e LICENSEID=... when running the container, or mount your own gurobi.lic
# at /gurobi.lic.
#
# The upstream repo (pmccaffrey6/immunoinformatics_platform) inherited
# committed credentials; this fork removes them.
ENV PATH="${PATH}:/gurobi952/linux64/bin"
ENV GRB_LICENSE_FILE=/gurobi.lic
ENV GUROBI_HOME=/gurobi952/linux64/

RUN wget https://packages.gurobi.com/9.5/gurobi9.5.2_linux64.tar.gz && \
    tar xzvf gurobi9.5.2_linux64.tar.gz && \
    conda create -y -n jessev --file /JessEV/packages.txt -c conda-forge
