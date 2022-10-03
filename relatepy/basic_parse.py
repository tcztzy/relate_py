import pandas as pd

def read_coal(filename):
    with open(filename) as f:
        groups = f.readline().split()
        epochs = [float(epoch) for epoch in f.readline().split()]
        coal = pd.read_table(f, sep=r"\s+", na_values="nan", names=["group1", "group2", *epochs])
    if len(groups) >= 1 + coal["group1"].max() and len(groups) >= 1 + coal["group2"].max():
        coal["group1"] = [groups[i] for i in coal["group1"]]
        coal["group2"] = [groups[i] for i in coal["group2"]]
    return coal.melt(id_vars=["group1", "group2"], var_name="epoch.start", value_name="haploid.coalescence.rate")
