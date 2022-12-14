#!/usr/bin/env python
from optparse import OptionParser
import pandas as pd

####################################################################################################
###     create cumulative region distribution file similar to mosdepth regio.dist.txt file       ###
####################################################################################################
if __name__ == '__main__':

    # read option
    optparser = OptionParser(usage="Usage:{} --bed coverage.per-base.bed".format(__file__))
    optparser.add_option('-b', '--bed', help='Per base coverage bed for panel regions', action='store')
    (options, args) = optparser.parse_args()
    if not options.bed:
        optparser.error('NO BED FILE GIVEN!')

    # create dataframe
    df = pd.read_csv(options.bed, sep="\t", names=["chr", "start", "stop", "depth"])

    # calculate length of region
    df['length']=df['stop'] - df['start']

    # sum length with equal depth
    df2 = pd.DataFrame(data=df.groupby('depth')['length'].sum())

    # sort by depth
    df2 = df2.sort_values(by='depth', ascending=False)

    # calculate cummulative percentage
    df2['cum_percent'] = df2.cumsum() / df2.sum()

    # print output + add missing values for every depth
    # 3 numbers after comma (= more specific than real mosdepth dist.txt file)
    for value in list(reversed(range(0, df2.index.tolist()[0] + 1, 1))):
        try:
            running_cum_percent = df2.loc[value, 'cum_percent']
            print(f'total\t{value}\t{running_cum_percent:.3f}')
        except KeyError:
            print(f'total\t{value}\t{running_cum_percent:.3f}')
