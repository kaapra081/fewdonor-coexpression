#!/usr/bin/env bash
# Download the two public CELLxGENE datasets used in the paper and verify checksums.
set -euo pipefail
mkdir -p data && cd data
# SEA-AD MTG Microglia-PVM (CELLxGENE dataset c76098ba-eed3-45b1-98f2-96fcac55ed18; dataset version ab53c08f-...)
curl -L -o seaad_mtg_microglia.h5ad https://datasets.cellxgene.cziscience.com/ab53c08f-03a3-4e0e-973a-44beabd78f60.h5ad
# Lee et al. 2026 fresh ex vivo myeloid, human PFC (CELLxGENE dataset 440da0a0-f321-4f40-9080-bb04031681e5; version 21a0c6d5-...)
curl -L -o lee2026_fresh_myeloid.h5ad https://datasets.cellxgene.cziscience.com/21a0c6d5-658c-434e-a479-1e51504aec5c.h5ad
shasum -a 256 -c ../data/checksums.sha256
