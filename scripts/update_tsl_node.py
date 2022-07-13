#!/usr/bin/env python

import subprocess
import pandas as pd
import tempfile
from pyminimc import util, pdf
import xml.etree.ElementTree as ET

from os.path import join


# read coarsened DataFrames
alpha_beta_partitions = {
    "alpha_partitions": [
        pd.DataFrame(
            pd.read_hdf(
                f"/Users/atumulak/Developer/minimc/data/tsl/endfb8/dfs/alpha_{i}.hdf5"
            )
        )
        for i in range(4)
    ],
    "beta_partitions": [
        pd.DataFrame(
            pd.read_hdf(
                f"/Users/atumulak/Developer/minimc/data/tsl/endfb8/dfs/beta_{i}.hdf5"
            )
        )
        for i in range(4)
    ],
}


# load input file to be modified
inputfilepath = "/Users/atumulak/Developer/minimc/benchmarks/broomstick.xml"
tree = ET.parse(inputfilepath)
tsl_node_path = "nuclides/continuous/nuclide/neutron/scatter/tsl"
if (tsl_node := tree.find(tsl_node_path)) is None:
    raise RuntimeError(f"node not found: {tsl_node_path}")

with tempfile.TemporaryDirectory() as tmpdir:
    print(f"Working in temporary directory: {tmpdir}")
    # modify input file
    for name, V_attribute in zip(["alpha", "beta"], ["beta_T", "E_T"]):
        partitions_node_name = f"{name}_partitions"
        if (partitions_node := tsl_node.find(partitions_node_name)) is not None:
            # remove existing partitions
            tsl_node.remove(partitions_node)
        # add empty partitions node
        partitions_node = ET.Element(partitions_node_name)
        # factorize full_df and save as separate DataFrames
        for i, full_df in enumerate(
            alpha_beta_partitions[partitions_node_name]
        ):
            partition_node = ET.Element("partition")
            for attribute, df in zip(
                ["CDF", "S", V_attribute], util.to_svd_dfs(full_df)
            ):
                df_path = join(tmpdir, f"{name}_{i}_{attribute}.hdf5")
                # save to temporary directory
                df.to_hdf(df_path, "pandas")
                # add path to DataFrame in partition node
                partition_node.set(attribute, df_path)
            partitions_node.append(partition_node)
        tsl_node.append(partitions_node)
    with tempfile.NamedTemporaryFile(dir=tmpdir) as modified_inputfile:
        tree.write(modified_inputfile)
        modified_inputfile.flush()  # stackoverflow.com/a/9422590/5101335
        # run minimc and parse output
        subprocess.run(
            [
                "/Users/atumulak/Developer/minimc/build/src/runminimc",
                modified_inputfile.name,
            ]
        )
        energy_node_path = "problemtype/fixedsource/energy/constant"
        if (energy_node := tree.find(energy_node_path)) is None or (
            energy := energy_node.get("energy")
        ) is None:
            raise RuntimeError(f"missing node or attribute: {energy_node_path}")
        temperature_node_path = "temperature/constant"
        if (temperature_node := tree.find(temperature_node_path)) is None or (
            temperature := temperature_node.get("c")
        ) is None:
            raise RuntimeError(
                f"missing node or attribute: {temperature_node_path}"
            )
        minimc_df = pdf.get_pdf_minimc(
            join(tmpdir, "minimc.out"), 1e6 * float(energy), float(temperature)
        )
        series = util.marginalize(minimc_df)
