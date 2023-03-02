import numpy as np
import pandas as pd
from collections import OrderedDict


def get_levels_dictionary(levels):
    temp_df = pd.read_csv(levels, sep="\t")
    uniqueLevels = temp_df["level"].unique()
    levelsDict = {elem: pd.DataFrame() for elem in uniqueLevels}

    for key in levelsDict.keys():
        levelsDict[key] = temp_df[:][temp_df["level"] == key]

    lvl1 = (
        levelsDict[1]["POS"].values,
        levelsDict[1]["REF"].values,
        levelsDict[1]["ALT"].values,
        levelsDict[1]["lineage"].values,
    )
    lvl2 = (
        levelsDict[2]["POS"].values,
        levelsDict[2]["REF"].values,
        levelsDict[2]["ALT"].values,
        levelsDict[2]["lineage"].values,
    )
    lvl3 = (
        levelsDict[3]["POS"].values,
        levelsDict[3]["REF"].values,
        levelsDict[3]["ALT"].values,
        levelsDict[3]["lineage"].values,
    )
    lvl4 = (
        levelsDict[4]["POS"].values,
        levelsDict[4]["REF"].values,
        levelsDict[4]["ALT"].values,
        levelsDict[4]["lineage"].values,
    )
    lvl5 = (
        levelsDict[5]["POS"].values,
        levelsDict[5]["REF"].values,
        levelsDict[5]["ALT"].values,
        levelsDict[5]["lineage"].values,
    )
    lineage4_pos = np.concatenate(
        (
            levelsDict[1].loc[levelsDict[1]["lineage"] == "L4", "POS"].values,
            levelsDict[2].loc[levelsDict[2]["lineage"] == "L4.9", "POS"].values,
        )
    )
    pos_all = np.setdiff1d(
        np.concatenate((lvl1[0], lvl2[0], lvl3[0], lvl4[0], lvl5[0])), lineage4_pos
    )

    return lvl1, lvl2, lvl3, lvl4, lvl5, lineage4_pos, pos_all


def get_variants_dataframe(variants_file, levels_file):
    chunks = pd.read_csv(variants_file, sep="\t", chunksize=50000, iterator=True)
    lineage4_pos, pos_all = get_levels_dictionary(levels_file)[5:7]
    chunk_list = []
    for chunk in chunks:
        result = chunk.melt(
            chunk[chunk.columns[chunk.columns.isin(["POS", "REF"])]],
            var_name="Sample",
            value_name="ALT",
        )
        result["Sample"] = result["Sample"].map(lambda x: str(x)[:-3])
        lineage4 = result.loc[result["POS"].isin(lineage4_pos)]
        lineage4 = lineage4[
            ~lineage4.duplicated(["Sample", "POS"], keep=False)
            | lineage4["ALT"].ne(lineage4["REF"])
        ]
        result = result.loc[result["POS"].isin(pos_all)]
        result = result.drop(result[result["REF"] == result["ALT"]].index)
        chunk_list.append(lineage4)
        chunk_list.append(result)
    df = pd.concat(chunk_list).reset_index(drop=True)
    return df


def count_level1_variants(myList):
    d = OrderedDict()
    for item in myList:
        if not len(item) == 0:
            caseless = item.casefold()
            try:
                d[caseless][1] += 1
            except KeyError:
                d[caseless] = [item, 1]

    myList = []
    for item, count in d.values():
        if not item.startswith("L8"):
            if count > 1:
                item = "{}".format(item)
            elif count == 1:
                item = "{} [{}]".format(item, "warning! only 1/2 snp is present")
        myList.append(item)
    return myList


def count_level2_variants(myList):
    d = OrderedDict()
    for item in myList:
        if not len(item) == 0:
            caseless = item.casefold()
            try:
                d[caseless][1] += 1
            except KeyError:
                d[caseless] = [item, 1]

    myList = []
    for item, count in d.values():
        if not item.startswith(("L2.2 (modern)", "L2.2 (ancient)")):
            if count > 1:
                item = "{}".format(item)
            elif count == 1:
                item = "{} [{}]".format(item, "warning! only 1/2 snp is present")
        elif item.startswith(("L2.2 (modern)", "L2.2 (ancient)")):
            if count > 1:
                item = "{}".format(item)
            elif count == 1:
                if not item.startswith("L2.2 (modern)"):
                    item = "{} [{}]".format(item, "warning! only 1/2 snp is present")
                if item.startswith("L2.2 (modern)"):
                    item = "{}".format(item)
        myList.append(item)
    return myList


def lineage2_decision(myList):
    lin2 = ["L2.2 (modern)", "L2.2 (ancient)"]
    altList = []

    for item in myList:
        if all(i in item for i in lin2) == True:
            item = list(set(item) - set(lin2))
            item.append("L2.2 (modern)")
        else:
            item = item

        altList.append(item)

    return altList


def lineage4_decision(call_list):
    lin4 = ["L4"]
    altList = []

    for item in call_list:
        if any(i in item for i in lin4):
            item = [x for x in item if x not in lin4]

        else:
            item.extend(["L4" for i in range(2)])

        altList.append(item)

    return altList


def lineage4_9_decision(call_list):
    lin4_9 = ["L4.9"]
    altList = []

    for item in call_list:
        if any(i in item for i in lin4_9):
            item = [x for x in item if x not in lin4_9]

        else:
            item.extend(["L4.9" for i in range(2)])

        altList.append(item)

    return altList


def barcoding(input_file, levels_file):
    lvl1, lvl2, lvl3, lvl4, lvl5 = get_levels_dictionary(levels_file)[0:5]
    df = get_variants_dataframe(input_file, levels_file)
    exp = (
        df["POS"].values[:, None],
        df["REF"].values[:, None],
        df["ALT"].values[:, None],
    )

    df["level_1"] = np.dot(
        np.logical_and.reduce(
            [
                np.equal(exp[0], lvl1[0]),
                np.equal(exp[1], lvl1[1]),
                np.equal(exp[2], lvl1[2]),
            ]
        ),
        lvl1[3],
    )

    df["level_2"] = np.dot(
        np.logical_and.reduce(
            [
                np.equal(exp[0], lvl2[0]),
                np.equal(exp[1], lvl2[1]),
                np.equal(exp[2], lvl2[2]),
            ]
        ),
        lvl2[3],
    )

    df["level_3"] = np.dot(
        np.logical_and.reduce(
            [
                np.equal(exp[0], lvl3[0]),
                np.equal(exp[1], lvl3[1]),
                np.equal(exp[2], lvl3[2]),
            ]
        ),
        lvl3[3],
    )

    df["level_4"] = np.dot(
        np.logical_and.reduce(
            [
                np.equal(exp[0], lvl4[0]),
                np.equal(exp[1], lvl4[1]),
                np.equal(exp[2], lvl4[2]),
            ]
        ),
        lvl4[3],
    )

    df["level_5"] = np.dot(
        np.logical_and.reduce(
            [
                np.equal(exp[0], lvl5[0]),
                np.equal(exp[1], lvl5[1]),
                np.equal(exp[2], lvl5[2]),
            ]
        ),
        lvl5[3],
    )

    df = df.drop(["REF", "ALT", "POS"], axis=1)
    df = df.replace("", np.nan)

    df = (
        df.groupby(["Sample"])
        .agg(lambda x: ",".join(x.dropna()))
        .reset_index()
        .reindex(columns=df.columns)
    )

    df2 = df.copy()
    df2["level_1"] = df2["level_1"].str.split(",")
    df2["level_2"] = df2["level_2"].str.split(",")
    df2["level_1"] = lineage4_decision(df2["level_1"])
    df2["level_2"] = lineage4_9_decision(df2["level_2"])
    df2["level_1"] = [count_level1_variants(item) for item in df2["level_1"]]
    df2["level_2"] = [count_level2_variants(item) for item in df2["level_2"]]
    df2["level_2"] = lineage2_decision(df2["level_2"])
    df2[["level_1", "level_2"]] = df2[["level_1", "level_2"]].applymap(
        lambda x: ", ".join(map(str, x))
    )
    df2 = df2.sort_values(
        by=["level_1", "level_2", "level_3", "level_4", "level_5"]
    ).reset_index(drop=True)

    return df2


result = barcoding(snakemake.input[0], snakemake.input[1])
result.to_csv(snakemake.output[0], sep="\t", index=False)
