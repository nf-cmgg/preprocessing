import argparse
import pybedtools
import pandas as pd
from pathlib import Path
import uuid
import csv


def main(per_base_bed: Path, panel_bed: Path) -> None:
    # create temporary intersection file
    panel_coverage_tmp = f"{uuid.uuid4()}.tmp"
    pybedtools.BedTool(per_base_bed).intersect(pybedtools.BedTool(panel_bed)).moveto(panel_coverage_tmp)

    # create dataframe
    df = pd.read_csv(panel_coverage_tmp, sep="\t", lineterminator="\n", names=["chr", "start", "stop", "depth"])

    # remove temporary file
    Path(panel_coverage_tmp).unlink()

    # calculate length of region
    df['length'] = df['stop'] - df['start']

    # sum length with equal depth
    df2 = pd.DataFrame(data=df.groupby('depth')['length'].sum())

    # sort by depth
    df2 = df2.sort_values(by='depth', ascending=False)

    # calculate cummulative percentage
    df2['cum_percent'] = df2.cumsum() / df2.sum()

    # print output + add missing values for every depth
    # 3 decimals (= more specific than real mosdepth dist.txt file)
    cum_region_dist = open(f'${meta.id}_${genelist_name}.region.dist.txt', 'w')
    with cum_region_dist as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for value in list(reversed(range(0, df2.index.tolist()[0]+1,1))):
            try:
                running_cum_percent = df2.loc[value, 'cum_percent']
                writer.writerow(["total", value, f"{running_cum_percent:.3f}"])
            except KeyError:
                writer.writerow(["total", value, f"{running_cum_percent:.3f}"])


if __name__ == '__main__':
    main(per_base_bed="$per_base_bed", panel_bed="$panel_bed")
