import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle

import pysam

def read_pair_generator(bam, refname, min_site, max_site):
    is_paired = {}
    for read in bam.fetch(refname, min_site, max_site+1):
        if read.query_name[-2] == '.':
            qname = read.query_name[:-1]
        else:
            qname = read.query_name

        if qname not in is_paired:
            is_paired[qname] = False
        else:
            is_paired[qname] = True

    read_dict = {}
    for read in bam.fetch(refname, min_site, max_site+1):
        if read.query_name[-2] == '.':
            qname = read.query_name[:-1]
        else:
            qname = read.query_name

        if not is_paired[qname]:
            yield read, None
            continue

        if qname not in read_dict:
            read_dict[qname] = read
        else:
            yield read_dict[qname], read
            del read_dict[qname]


def extract(query_mutations, input_bam, output, same_read):

    # Parse SNPs and Indels from query_mutations
    with open(query_mutations) as infile:
        lines = infile.read().splitlines()
        snps = []
        insertions = []
        deletions = []
        for line in lines:
            line = [s.strip() for s in line.split(',')]
            if any(n in line[0] for n in 'ACGT'):
                if ':' in line[0]:
                    insertions = line
                elif len(line) > 0:
                    snps = line
            elif ':' in line[0]:
                deletions = line
    try:
        if len(insertions) > 0:
            insertions = [
                (int(s.split(':')[0][1:]), s.split(':')[1][1:-2].strip('\''))
                for s in insertions
            ]
        if len(deletions) > 0:
            deletions = [
                (int(s.split(':')[0][1:]), int(s.split(':')[1][:-1]))
                for s in deletions
            ]
        snp_sites = [int(m[1:len(m)-1])-1 for m in snps if m]
        indel_sites = [s[0] for s in insertions if s] +\
            [s[0] for s in deletions if s]
    except Exception:
        print('extract: Error parsing', query_mutations)
        print('extract: See README for formatting requirements.')
        return -1

    all_sites = snp_sites + indel_sites
    all_sites.sort()
    snp_dict = {int(mut[1:len(mut)-1])-1: mut[-1] for mut in snps if mut}

    # Open input_bam and output for reading/writing
    print("extract: Extracting read pairs with specified mutations")
    try:
        samfile = pysam.AlignmentFile(input_bam, 'rb')
    except ValueError:
        print('extract: Missing index file, try running samtools index',
              input_bam)
        return -1
    refname = samfile.get_reference_name(0)

    outfile = pysam.AlignmentFile(output, 'wb', template=samfile)

    min_site = min(all_sites) - 150
    if min_site < 0:
        min_site = 0
    max_site = max(all_sites) + 150
    if max_site > samfile.get_reference_length(refname):
        max_site = samfile.get_reference_length(refname)

    reads_considered = []

    for read1, read2 in read_pair_generator(samfile, refname, min_site,
                                            max_site):
        indels_found = {mut: False for mut in insertions + deletions}
        snps_found = {site: False for site in snp_dict}

        for x in [read1, read2]:

            if x is None:
                continue

            if any([val is None for val in [x.query_alignment_sequence,
                                            x.query_alignment_qualities,
                                            x.cigarstring]]):
                continue

            ref_pos = set(x.get_reference_positions())
            start = x.reference_start
            sites_in = list(ref_pos & set(snp_sites))

            seq = x.query_alignment_sequence

            cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)

            # Find insertions
            if 'I' in x.cigarstring:
                i = 0
                del_spacing = 0
                ins_spacing = 0
                ranges = []

                for m in cigar:
                    if m[1] == 'I':
                        current_ins = (start+i+del_spacing-ins_spacing,
                                       seq[i:i+int(m[0])])
                        if current_ins in indels_found:
                            indels_found[current_ins] = True
                        i += int(m[0])
                        ins_spacing += int(m[0])
                    elif m[1] == 'D':
                        del_spacing += int(m[0])
                    elif m[1] == 'M':
                        ranges.append([i, i+int(m[0])])
                        i += int(m[0])

                seq = ''.join([seq[r[0]:r[1]] for r in ranges])

            # Find deletions
            if 'D' in x.cigarstring:
                i = x.reference_start
                for m in cigar:
                    if m[1] == 'M':
                        i += int(m[0])
                    elif m[1] == 'D':
                        current_del = (i, int(m[0]))
                        if current_del in deletions:
                            indels_found[current_del] = True
                        i += int(m[0])

            # Find SNPs
            c_dict = dict(zip(x.get_reference_positions(), seq))
            for s in sites_in:
                if snp_dict[s] == c_dict[s]:
                    snps_found[s] = True

        if same_read:  # check that ALL muts are present
            mut_missing = False
            for mut in indels_found:
                if not indels_found[mut]:
                    mut_missing = True
                    break

            for site in snps_found:
                if not snps_found[site]:
                    mut_missing = True
                    break

            if not mut_missing:
                if read1 not in reads_considered:
                    outfile.write(read1)
                    reads_considered.append(read1)
                if read2 is not None:
                    outfile.write(read2)
                    reads_considered.append(read2)

        else:  # at least 1 mutation present
            mut_found = False
            for mut in indels_found:
                if indels_found[mut]:
                    mut_found = True
                    break

            for site in snps_found:
                if snps_found[site]:
                    mut_found = True
                    break

            if mut_found:
                if read1 not in reads_considered:
                    outfile.write(read1)
                    reads_considered.append(read1)
                if read2 is not None:
                    if read2 not in reads_considered:
                        outfile.write(read2)
                        reads_considered.append(read2)

    samfile.close()
    outfile.close()

    print(f'extract: Output saved to {output}')
    return reads_considered


def filter(query_mutations, input_bam, min_site, max_site, output):

    # Load data
    with open(query_mutations) as infile:
        lines = infile.read().splitlines()
        snps = []
        insertions = []
        deletions = []
        for line in lines:
            line = [s.strip() for s in line.split(',')]
            if any(n in line[0] for n in 'ACGT'):
                if ':' in line[0]:
                    insertions = line
                elif len(line) > 0:
                    snps = line
            elif ':' in line[0]:
                deletions = line
    try:
        # parse tuples from indels
        if len(insertions) > 0:
            insertions = [(int(s.split(':')[0][1:]),
                           s.split(':')[1][1:-2].strip('\''))
                          for s in insertions]
        if len(deletions) > 0:
            deletions = [(int(s.split(':')[0][1:]), int(s.split(':')[1][:-1]))
                         for s in deletions]

        # get loci for all mutations
        snp_sites = [int(m[1:len(m)-1])-1 for m in snps if m]
        indel_sites = [s[0] for s in insertions if s] +\
            [s[0] for s in deletions if s]
    except ValueError:
        print('filter: Error parsing', query_mutations)
        print('filter: See README for formatting requirements.')
        return -1

    snp_dict = {int(mut[1:len(mut)-1])-1: mut[-1] for mut in snps if mut}
    reads_considered = []

    try:
        samfile = pysam.AlignmentFile(input_bam, 'rb')
    except ValueError:
        print('filter: Missing index file, try running samtools index',
              input_bam)
        return -1
    print("filter: Filtering out reads with specified mutations")
    refname = samfile.get_reference_name(0)
    all_sites = snp_sites + indel_sites
    all_sites.sort()

    for site in all_sites:
        itr = samfile.fetch(refname, site, site+1)

        for x in itr:
            if x is None:
                continue

            if any([val is None for val in [x.query_alignment_sequence,
                                            x.query_alignment_qualities,
                                            x.cigarstring]]):
                continue

            ref_pos = set(x.get_reference_positions())
            start = x.reference_start
            sites_in = list(ref_pos & set(snp_sites))

            seq = x.query_alignment_sequence

            cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)

            # note: inserts add to read length, but not alignment coordinate
            if 'I' in x.cigarstring:
                i = 0
                del_spacing = 0
                ins_spacing = 0
                ranges = []
                insert_found = False
                for m in cigar:
                    if m[1] == 'I':
                        if (start+i+del_spacing-ins_spacing,
                                seq[i:i+int(m[0])]) in insertions:
                            reads_considered.append(x.query_name)
                            insert_found = True
                            break
                        i += int(m[0])
                        ins_spacing += int(m[0])
                    elif m[1] == 'D':
                        del_spacing += int(m[0])
                        continue
                    elif m[1] == 'M':
                        ranges.append([i, i+int(m[0])])
                        i += int(m[0])

                seq = ''.join([seq[r[0]:r[1]] for r in ranges])

                if insert_found:
                    continue

            if 'D' in x.cigarstring:
                i = x.reference_start
                for m in cigar:
                    if m[1] == 'M':
                        i += int(m[0])
                    elif m[1] == 'D':
                        if (i, int(m[0])) in deletions:
                            reads_considered.append(x.query_name)
                            continue
                        i += int(m[0])
            c_dict = dict(zip(x.get_reference_positions(), seq))
            for s in sites_in:
                if snp_dict[s] == c_dict[s]:
                    reads_considered.append(x.query_name)
                    break
    samfile.close()

    # Run again, this time also getting the paired read
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    outfile = pysam.AlignmentFile(output, 'wb', template=samfile)
    final_reads = []

    itr = samfile.fetch(refname, int(min_site), int(max_site)+1)

    for x in itr:
        if x.query_name in reads_considered:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)

    samfile.close()
    outfile.close()

    print(f'filter: Output saved to {output}')

    return final_reads


def plot_covariants(covariants, output, num_clusters,
                    min_mutations, nt_muts):

    covariants_df = pd.read_csv(covariants, sep='\t')

    covariants_df['nt_mutations'] = covariants_df['nt_mutations'].apply(
        lambda x: x if len(x.split(' ')) >= min_mutations else pd.NA
    )
    covariants_df = covariants_df.dropna(subset=['nt_mutations'])

    nt_to_aa = {}
    for idx, row in covariants_df.iterrows():
        for i, mut in enumerate(row['nt_mutations'].split(' ')):
            if mut not in nt_to_aa:
                if row['aa_mutations'].split(' ')[i] != 'NA':
                    nt_to_aa[mut] = row['aa_mutations'].split(' ')[i]
                else:
                    nt_to_aa[mut] = mut

    # Merge rows with identical covariants and add their counts
    covariants_df = covariants_df \
        .groupby(covariants_df['nt_mutations'].astype(str)) \
        .agg(
            {
                'cluster_depth': 'sum',
                'coverage_start': 'max',
                'coverage_end': 'min',
                'frequency': 'sum'
            }
        ) \
        .reset_index() \
        .sort_values('frequency', ascending=False)

    covariants_df['Coverage_ranges'] = list(
        zip(
            covariants_df['coverage_start'],
            covariants_df['coverage_end']
        )
    )
    covariants_df = covariants_df.head(num_clusters)
    coverage_ranges = covariants_df['Coverage_ranges'].tolist()
    log_frequency = np.log10(covariants_df['frequency']).tolist()

    # Get all unique mutations found in the sample
    unique_mutations = set()
    for sublist in covariants_df['nt_mutations']:
        unique_mutations.update(sublist.split(' '))

    def get_position(mut, nt_muts):
        # if nt_muts:
        if '-' in mut:
            return int(mut[1:-1].split('-')[0])
        if '+' in mut:
            return int(mut[1:-1].split('+')[0])
        m = re.search(r'(\d+)', mut)
        if m:
            return int(m.group(1))
        return 0


    # Populate dataframe with cluster frequencies
    sites = {mut: get_position(mut, nt_muts) for mut in unique_mutations}
    colnames = list(unique_mutations)
    colnames.sort(key=lambda x: sites[x])

    index = [f'{i}' for i in range(len(covariants_df))]
    plot_df = pd.DataFrame(columns=colnames, index=index)
    for idx, cov in enumerate(covariants_df['nt_mutations']):
        sample_row = np.empty(len(colnames))
        sample_row[:] = np.nan
        coverage_range = coverage_ranges[idx]
        muts = [m for m in str(cov).split() if m]
        for col in colnames:
            if sites[col] in range(coverage_range[0], coverage_range[1]):
                # Placeholder value for reference bases
                sample_row[colnames.index(col)] = 1.0
            for mut in muts:
                # match whole mutation tokens, not characters
                if col == mut:
                    sample_row[colnames.index(col)] = log_frequency[idx]
                    break
        plot_df.loc[index[idx]] = sample_row

    plot_df = plot_df.astype(float)

    if not nt_muts:
        plot_df.columns = [nt_to_aa[col] for col in plot_df.columns]
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(15, 10))

    colors = ['#FFFFB2', '#FECC5C', '#FD8D3C', '#E31A1C']
    cmap = mcolors.LinearSegmentedColormap.from_list(name='custom',
                                                     colors=colors)
    ax = sns.heatmap(plot_df, cmap=cmap,
                     cbar_kws={'label': 'log10 Frequency', 'shrink': 0.5},
                     vmax=0,
                     linewidths=1.5, linecolor='white', square=True)

    gray = mcolors.LinearSegmentedColormap.from_list(
        name='custom', colors=['#BEBEBE', '#BEBEBE'])
    plot = sns.heatmap(plot_df, cmap=gray, vmin=0, vmax=2,
                       linewidths=1.75, linecolor='white', square=True,
                       mask=plot_df <= 0, cbar=False, ax=ax)

    # Add line between adjacent non-NaN cells
    for i, row in enumerate(plot_df.values):
        j = 0
        while j < len(row):
            if not np.isnan(row[j]):

                if j > 0 and not np.isnan(row[j - 1]):
                    line = plt.Line2D((j, j), (i, i + 1),
                                      color='gray', linewidth=2.75, zorder=1)
                    ax.add_artist(line)
            j += 1

    # Add rectangles for non-NaN cells
    for i, row in enumerate(plot_df.values):
        for j, val in enumerate(row):
            if val <= 1:
                width = np.count_nonzero(~np.isnan(row))
                rect = Rectangle((j, i), width, 1, fill=False,
                                 edgecolor='black', linewidth=4)
                ax.add_patch(rect)
                break

    for _, spine in plot.spines.items():
        spine.set_visible(True)
    plt.yticks([])

    legend_elements = [
        Patch(facecolor='#BEBEBE', edgecolor='black',
              label='Reference sequence'),
        Patch(facecolor='white', edgecolor='black', label='No coverage'),
    ]
    ax.legend(handles=legend_elements,
              bbox_to_anchor=(1.20, 0.40), loc="lower left", borderaxespad=0)

    fig.tight_layout()
    plt.savefig(output, bbox_inches='tight')

    print(f'plot-covariants: Output saved to {output}')
