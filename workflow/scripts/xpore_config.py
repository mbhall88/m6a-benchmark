import sys
sys.stderr = open(snakemake.log, "w")

from collections import defaultdict
from pathlib import Path

import yaml

data = defaultdict()
for p in map(Path, snakemake.input.jsons):
    cond, rep = p.parent.name.split("_")
    data[cond][rep] = p.resolve()

criteria = {
    "readcount_min": snakemake.params.readcount_min,
    "readcount_max": snakemake.params.readcount_max,
}

yml = {"out": snakemake.params.outdir, "data": data, "criteria": criteria}

with open(snakemake.output.configuration, "w") as ofp:
    yaml.dump(yml, ofp, default_flow_style=False)
