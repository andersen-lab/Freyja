import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle
import concurrent.futures as cf
from concurrent.futures import ProcessPoolExecutor
import pysam
from Bio.Seq import MutableSeq
from Bio import SeqIO

from freyja.read_analysis_utils import nt_position, get_colnames_and_sites, \
    read_pair_generator, \
    filter_covariants_output, parse_gff, get_gene


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


def run_parallel(cores, regions, input_bam, ref_fasta, gff_file, min_quality,
                 min_count, spans_region):
    with ProcessPoolExecutor(max_workers=cores) as p:
        futures = {p.submit(process_covariants, input_bam,
                            region[0], region[1], ref_fasta, gff_file,
                            min_quality, min_count, spans_region)
                   for region in regions}
        for future in cf.as_completed(futures):
            yield future.result()


def covariants(input_bam, min_site, max_site, output,
               ref_fasta, gff_file, min_quality, min_count, spans_region,
               sort_by, cores):
    # Split regions into smaller chunks for parallel processing
    regions = [(min_site, max_site)]
    if cores > 1:
        region_size = (max_site - min_site) // cores
        regions = [(min_site + i * region_size, min_site +
                    (i + 1) * region_size) for i in range(cores)]
        regions[-1] = (regions[-1][0], max_site)

    # Run parallel processing
    results = run_parallel(cores, regions, input_bam, ref_fasta,
                           gff_file, min_quality, min_count, spans_region)

    # Aggregate results

    results = list(results)
    df = pd.concat(results, ignore_index=True)

    df = df.sort_values('Max_count', ascending=False).drop_duplicates(
        subset='Covariants', keep='first')

    # Sort patterns
    if sort_by.lower() == 'count':
        df = df.sort_values('Count', ascending=False)
    elif sort_by.lower() == 'freq':
        df = df.sort_values('Freq', ascending=False)
    elif sort_by.lower() == 'site':
        df['sort_col'] = [nt_position(s.split(' ')[0]) for s in df.Covariants]
        df = df.sort_values('sort_col').drop(labels='sort_col', axis=1)

    df.to_csv(output, sep='\t', index=False)
    print(f'covariants: Output saved to {output}')
    return df


def process_covariants(input_bam, min_site, max_site, ref_fasta, gff_file,
                       min_quality, min_count, spans_region):

    # Load reference genome
    ref_genome = MutableSeq(next(SeqIO.parse(ref_fasta, 'fasta')).seq)

    # Load gene annotations (if present)
    if gff_file is not None:
        gene_positions = parse_gff(gff_file)

        # Split ORF1ab for SARS-CoV-2
        if 'ORF1ab' in gene_positions:
            del gene_positions['ORF1ab']
            gene_positions['ORF1a'] = (266, 13468)
            gene_positions['ORF1b'] = (13468, 21555)
        elif 'orf1ab' in gene_positions:
            del gene_positions['orf1ab']
            gene_positions['orf1a'] = (266, 13468)
            gene_positions['orf1b'] = (13468, 21555)

    # Open input bam file for reading
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    refname = samfile.get_reference_name(0)

    co_muts = {}
    co_muts_by_gene = {}
    coverage = {}
    co_muts_region = {}

    # Check if index file exists
    try:
        samfile.fetch(refname, min_site, max_site+1)
    except ValueError:
        print((f'covariants: Input bamfile missing corresponding index file. '
               f'Try running:\n samtools index {input_bam}'))
        return -1

    for read1, read2 in read_pair_generator(samfile, refname, min_site,
                                            max_site+1):

        coverage_start = None
        coverage_end = None
        co_mut_start = None
        co_mut_end = None
        nt_to_aa = {}

        muts_final = []
        for x in [read1, read2]:
            if x is None:
                continue

            if any([val is None for val in [x.query_alignment_sequence,
                                            x.query_alignment_qualities,
                                            x.cigarstring]]):
                continue

            start = x.reference_start
            end = x.reference_end
            seq = x.query_alignment_sequence
            quals = x.query_alignment_qualities
            seq = ''.join(
                [seq[i] if quals[i] >= min_quality else 'N'
                 for i in range(len(seq))])

            insertions_found = []
            ins_offsets = {}
            ins_offset = 0
            last_ins_site = 0

            deletions_found = []
            del_offsets = {}
            del_offset = 0
            last_del_site = 0

            snps_found = []

            # Update coverage ranges
            if coverage_start is None or start < coverage_start:
                coverage_start = start

            if coverage_end is None or end > coverage_end:
                coverage_end = end

            cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)
            if 'I' in x.cigarstring:
                i = 0
                del_spacing = 0
                ins_spacing = 0

                for m in cigar:
                    if m[1] == 'I':
                        insertion_site = start+i+del_spacing-ins_spacing
                        insertion_len = int(m[0])

                        ins_offsets[(last_ins_site, insertion_site)
                                    ] = ins_offset
                        last_ins_site = insertion_site
                        ins_offset += insertion_len
                        if 'N' in seq[i:i+insertion_len]:
                            continue

                        insertions_found.append(
                            (insertion_site,
                             seq[i:i+insertion_len])
                        )
                        # only insertion site in reference genome
                        if co_mut_start is None \
                                or insertion_site < co_mut_start:
                            co_mut_start = insertion_site
                        if co_mut_end is None \
                                or insertion_site > co_mut_end:
                            co_mut_end = insertion_site

                        i += insertion_len
                        ins_spacing += insertion_len
                    elif m[1] == 'D':
                        del_spacing += int(m[0])
                    elif m[1] == 'M':
                        i += int(m[0])

                current_ins_site = start+i+del_spacing-ins_spacing
                ins_offsets[(last_ins_site, current_ins_site)] = ins_offset
                last_ins_site = current_ins_site

            if 'D' in x.cigarstring:
                i = 0
                for m in cigar:
                    if m[1] == 'M':
                        i += int(m[0])
                    elif m[1] == 'D':
                        deletion_site = start+i
                        deletion_len = int(m[0])

                        deletions_found.append((deletion_site, deletion_len))
                        if co_mut_start is None \
                                or deletion_site < co_mut_start:
                            co_mut_start = deletion_site
                        if co_mut_end is None \
                                or deletion_site + deletion_len > co_mut_end:
                            co_mut_end = deletion_site + deletion_len

                        del_offsets[(last_del_site, start+i)] = del_offset
                        last_del_site = deletion_site
                        del_offset += deletion_len

                        i += deletion_len
                del_offsets[(last_del_site, start+i)] = del_offset
                last_del_site = start+i

            # Find SNPs

            softclip_offset = 0
            if cigar[0][1] == 'S':
                softclip_offset += int(cigar[0][0])
            if (len(cigar) > 1 and cigar[0][1] == 'H' and cigar[1][1] == 'S'):
                softclip_offset += int(cigar[1][0])

            pairs = x.get_aligned_pairs(matches_only=True)

            for tup in pairs:
                read_site, ref_site = tup

                read_site -= softclip_offset
                ref_base = ref_genome[ref_site]

                if seq[read_site] != 'N' and ref_base != seq[read_site]:
                    snps_found.append(
                        f'{ref_base.upper()}{ref_site+1}{seq[read_site]}'
                    )
                    if co_mut_start is None or ref_site + 1 < co_mut_start:
                        co_mut_start = ref_site + 1
                    if co_mut_end is None or ref_site + 1 > co_mut_end:
                        co_mut_end = ref_site + 1

            # Get corresponding amino acid mutations
            if gff_file is not None:
                for ins in insertions_found:
                    locus = ins[0]
                    gene_info = get_gene(locus, gene_positions)
                    ins_string = str(ins).replace(' ', '')
                    if gene_info is None or ins_string in nt_to_aa:
                        continue

                    gene, start_site = gene_info
                    aa_locus = ((locus - start_site) // 3) + 2

                    if len(ins[1]) % 3 == 0:
                        insertion_seq = MutableSeq(ins[1]).translate()
                    else:
                        insertion_seq = ''

                    aa_mut = (f'{ins_string}({gene}:INS{aa_locus}'
                              f'{insertion_seq})')
                    nt_to_aa[ins_string] = aa_mut

                for deletion in deletions_found:
                    locus = deletion[0]
                    gene_info = get_gene(locus, gene_positions)
                    deletion_string = str(deletion).replace(' ', '')
                    if gene_info is None or deletion_string in nt_to_aa:
                        continue

                    gene, start_site = gene_info
                    aa_locus = ((locus - start_site) // 3) + 2

                    del_length = deletion[1] // 3
                    if del_length > 1:
                        aa_mut = (
                            f'{deletion_string}({gene}:DEL{aa_locus}/'
                            f'{aa_locus+del_length-1})'
                        )
                    else:
                        aa_mut = f'{deletion_string}({gene}:DEL{aa_locus})'
                    nt_to_aa[deletion_string] = aa_mut

                for snp in snps_found:
                    locus = int(snp[1:-1])
                    gene_info = get_gene(locus, gene_positions)
                    if gene_info is None or snp in nt_to_aa:
                        continue

                    gene, start_site = gene_info
                    codon_position = (locus - start_site) % 3
                    aa_locus = ((locus - codon_position - start_site) // 3) + 1

                    ref_codon = ref_genome[locus - codon_position - 1:
                                           locus - codon_position + 2]
                    ref_aa = ref_codon.translate()

                    # Adjust for indels
                    ins_offset = 0
                    for r in ins_offsets:
                        if locus in range(r[0], r[1]) or locus == r[1]:
                            ins_offset = ins_offsets[r]
                    del_offset = 0
                    for r in del_offsets:
                        if locus in range(r[0], r[1]) or locus == r[1]:
                            del_offset = del_offsets[r]

                    read_start = (locus - start) - codon_position + \
                        ins_offset - del_offset - 1
                    read_end = (locus - start) - codon_position + \
                        ins_offset - del_offset + 2

                    alt_codon = MutableSeq(seq[read_start:read_end])

                    if len(alt_codon) % 3 != 0 or len(alt_codon) == 0 \
                            or 'N' in alt_codon:
                        continue
                    alt_aa = alt_codon.translate()

                    aa_mut = f'{snp}({gene}:{ref_aa}{aa_locus}{alt_aa})'
                    nt_to_aa[snp] = aa_mut

            for mut in insertions_found + deletions_found + snps_found:
                muts_final.append(str(mut).replace(' ', ''))

        # Skip update if reads do not span region
        if spans_region:
            if coverage_start > min_site or coverage_end < max_site:
                continue

        muts_final = sorted(list(set(muts_final)), key=nt_position)
        if len(muts_final) == 0:
            continue

        # Update
        key = ' '.join(muts_final)
        if key not in co_muts:
            co_muts[key] = 1
            coverage[key] = [coverage_start, coverage_end]
            co_muts_region[key] = (co_mut_start, co_mut_end)
        else:
            co_muts[key] += 1

            # Update coverage to intersection of read pairs with covariants
            if coverage_start > coverage[key][0]:
                coverage[key][0] = coverage_start
            if coverage_end < coverage[key][1]:
                coverage[key][1] = coverage_end

        # Should only have to do this when key not in co_muts
        # but aa_mut can sometimes change for snps due to indels
        if gff_file is not None:
            muts_final_aa = []
            for mut in muts_final:
                if mut in nt_to_aa:
                    muts_final_aa.append(nt_to_aa[mut])
                else:
                    muts_final_aa.append(mut)

            aa_string = ' '.join(muts_final_aa)
            if key not in co_muts_by_gene \
                    or len(co_muts_by_gene[key]) < len(aa_string):
                co_muts_by_gene[key] = aa_string

    samfile.close()

    # Go through samfile another time to determine number of read pairs
    # for positions covered by each set of covariants
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    samfile.fetch(refname, min_site, max_site+1)
    co_muts_max_reads = {}

    for read1, read2 in read_pair_generator(samfile, refname, min_site,
                                            max_site+1):

        coverage_start = None
        coverage_end = None

        for x in [read1, read2]:

            if x is None:
                continue

            if any([val is None for val in [x.query_alignment_sequence,
                                            x.query_alignment_qualities,
                                            x.cigarstring]]):
                continue

            start = x.reference_start
            end = x.reference_end
            seq = x.query_alignment_sequence
            quals = x.query_alignment_qualities
            seq = ''.join(
                [seq[i] if quals[i] >= min_quality else 'N'
                 for i in range(len(seq))])

            # Update coverage ranges
            if coverage_start is None or start < coverage_start:
                coverage_start = start
            if coverage_end is None or end > coverage_end:
                coverage_end = end

        # Skip update if reads do not span region
        if spans_region:
            if coverage_start > min_site or coverage_end < max_site:
                continue

        if coverage_start is None or coverage_end is None:
            continue

        for key in co_muts_region:
            co_mut_start, co_mut_end = co_muts_region[key]

            if coverage_start <= co_mut_start and coverage_end >= co_mut_end:
                if key not in co_muts_max_reads:
                    co_muts_max_reads[key] = 1
                else:
                    co_muts_max_reads[key] += 1

    samfile.close()

    # Aggregate dictionaries into dataframe
    df = pd.DataFrame()
    if gff_file is not None:
        df['Covariants'] = [co_muts_by_gene[k] for k in co_muts]
    else:
        df['Covariants'] = [k for k in co_muts]
    df['Count'] = [co_muts[k] for k in co_muts]
    df['Max_count'] = [co_muts_max_reads[k] for k in co_muts]
    df['Freq'] = [co_muts[k] / co_muts_max_reads[k] for k in co_muts]
    df['Coverage_start'] = [coverage[k][0] for k in co_muts]
    df['Coverage_end'] = [coverage[k][1] for k in co_muts]

    df = df[df['Count'] >= min_count]

    return df


def plot_covariants(covariants, output, num_clusters,
                    min_mutations, nt_muts, vmin, vmax):

    covariants_df = pd.read_csv(covariants, sep='\t', header=0)

    if not any(':' in mut for mut in covariants_df['Covariants'][0]):
        nt_muts = True

    # Clean up covariants output
    covariants_df['Covariants'] = covariants_df['Covariants'] \
        .str.replace(', ', ',') \
        .str.split(' ') \
        .apply(filter_covariants_output, args=(nt_muts, min_mutations))

    # Drop rows with NA values and keep top n clusters
    covariants_df = covariants_df.dropna()

    # Merge rows with identical covariants and add their counts
    covariants_df = covariants_df \
        .groupby(covariants_df['Covariants'].astype(str)) \
        .agg(
            {
                'Count': 'sum',
                'Coverage_start': 'max',
                'Coverage_end': 'min',
                'Freq': 'sum'
            }
        ) \
        .reset_index() \
        .sort_values('Freq', ascending=False)

    # Cast back to list
    covariants_df['Covariants'] = covariants_df['Covariants'] \
        .str.strip('][') \
        .str.replace("'", "") \
        .str.split(pat=', ')

    covariants_df['Coverage_ranges'] = list(
        zip(
            covariants_df['Coverage_start'],
            covariants_df['Coverage_end']
        )
    )
    covariants_df = covariants_df.head(num_clusters)
    coverage_ranges = covariants_df['Coverage_ranges'].tolist()
    log_frequency = np.log10(covariants_df['Freq']).tolist()

    # Get all unique mutations found in the sample
    unique_mutations = set()
    for sublist in covariants_df['Covariants']:
        unique_mutations.update(sublist)

    # Populate dataframe with cluster frequencies
    colnames, sites = get_colnames_and_sites(unique_mutations, nt_muts)
    index = [f'CP{i}' for i in range(len(covariants_df))]
    plot_df = pd.DataFrame(columns=colnames, index=index)
    for pattern in enumerate(covariants_df['Covariants']):
        sample_row = np.empty(len(colnames))
        sample_row[:] = np.nan
        coverage_range = coverage_ranges[pattern[0]]
        for col in colnames:
            if sites[col] in range(coverage_range[0], coverage_range[1]):
                # Placeholder value for refererence bases
                sample_row[colnames.index(col)] = 1.0
            for mut in pattern[1]:
                if col in mut:
                    sample_row[colnames.index(
                        col)] = log_frequency[pattern[0]]
        plot_df.loc[index[pattern[0]]] = sample_row

    plot_df = plot_df.astype(float)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(15, 10))

    colors = ['#FFFFB2', '#FECC5C', '#FD8D3C', '#E31A1C']
    cmap = mcolors.LinearSegmentedColormap.from_list(name='custom',
                                                     colors=colors)
    ax = sns.heatmap(plot_df, cmap=cmap,
                     cbar_kws={'label': 'log10 Frequency', 'shrink': 0.5},
                     vmin=vmin,
                     vmax=vmax,
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
