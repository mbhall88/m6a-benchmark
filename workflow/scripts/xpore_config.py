import sys
sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
from pathlib import Path

import yaml

data = defaultdict(dict)
for p in map(Path, snakemake.input.jsons):
    cond, rep = p.parent.name.split("_")
    data[cond][rep] = str(p.resolve().parent)

criteria = {
    "readcount_min": snakemake.params.readcount_min,
    "readcount_max": snakemake.params.readcount_max,
}

yml = {"out": str(snakemake.params.outdir), "data": dict(data), "criteria": criteria}

with open(snakemake.output.configuration, "w") as ofp:
    yaml.dump(yml, ofp, default_flow_style=False)
