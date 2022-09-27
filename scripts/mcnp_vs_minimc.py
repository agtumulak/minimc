#!/usr/bin/env python

import pandas as pd
from matplotlib import pyplot as plt
from pyminimc import util, estimator


mcnp = (
    estimator.MCNP(
        "/Users/atumulak/Developer/mcnp6-runs/single_zone/single_zone.mctal"
    )
    .marginalize("energy")
    .plot(label="MCNP")
)

util.run_minimc_with_tsl(
    "/Users/atumulak/Developer/minimc/build/src/runminimc",
    "/Users/atumulak/Developer/minimc/benchmarks/single_zone.xml",
    [
        pd.DataFrame(
            pd.read_hdf(
                f"/Users/atumulak/Developer/minimc/data/tsl/endfb8/dfs/beta_{i}.hdf5"
            )
        )
        for i in range(4)
    ],
    [
        pd.DataFrame(
            pd.read_hdf(
                f"/Users/atumulak/Developer/minimc/data/tsl/endfb8/dfs/alpha_{i}.hdf5"
            )
        )
        for i in range(4)
    ],
).marginalize("energy").plot(label="MiniMC")

plt.legend()
plt.show()
