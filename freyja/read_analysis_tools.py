import os
import re
import pysam


def extract(query_mutations, input_bam, output, refname):

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
        # parse tuples from indel strings
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

        # get loci for all mutations
        snp_sites = [int(m[1:len(m)-1])-1 for m in snps if m]
        indel_sites = [s[0] for s in insertions if s] +\
            [s[0] for s in deletions if s]
    except Exception:
        print('extract: Error parsing', query_mutations)
        print('extract: See README for formatting requirements.')
        return -1

    snp_dict = {int(mut[1:len(mut)-1])-1: mut[-1] for mut in snps if mut}

    try:
        samfile = pysam.AlignmentFile(input_bam, 'rb')
    except ValueError:
        print('extract: Missing index file, try running samtools index',
              input_bam)
        return -1

    reads_considered = []

    all_sites = snp_sites + indel_sites
    all_sites.sort()

    for site in all_sites:
        itr = samfile.fetch(refname, site, site+1)

        for x in itr:
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
            else:
                c_dict = dict(zip(x.get_reference_positions(), seq))
                for s in sites_in:
                    if snp_dict[s] == c_dict[s]:
                        reads_considered.append(x.query_name)
                        break
    samfile.close()

    # Run again, this time also getting the paired read
    samfile = pysam.AlignmentFile(input_bam, 'rb')

    outfile_name = input_bam.split('/')[-1][:-4] + '_extracted.bam'

    outfile_path = os.path.join(output, outfile_name)
    outfile = pysam.AlignmentFile(outfile_path, 'wb', template=samfile)
    final_reads = []

    itr = samfile.fetch(refname, min(all_sites), max(all_sites))

    for x in itr:
        if x.query_name not in reads_considered:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)

    samfile.close()
    outfile.close()

    print(f'Output saved to {outfile_path}')

    return final_reads


def filter(query_mutations, input_bam, min_site, max_site, output, refname):

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
    except ValueError:
        print('filter: Error parsing', query_mutations)
        print('filter: See README for formatting requirements.')
        return -1

    snp_dict = {int(mut[1:len(mut)-1])-1: mut[-1] for mut in snps if mut}
    reads_considered = []

    try:
        samfile = pysam.AlignmentFile(input_bam, 'rb')
    except ValueError:
        print('extract: Missing index file, try running samtools index',
              input_bam)
        return -1

    itr = samfile.fetch(refname, int(min_site), int(max_site)+1)

    for x in itr:
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
        else:
            c_dict = dict(zip(x.get_reference_positions(), seq))
            for s in sites_in:
                if snp_dict[s] == c_dict[s]:
                    reads_considered.append(x.query_name)
                    break
    samfile.close()

    # Run again, this time also getting the paired read
    samfile = pysam.AlignmentFile(input_bam, 'rb')

    outfile_name = input_bam.split('/')[-1][:-4] + '_filtered.bam'

    outfile_path = os.path.join(output, outfile_name)
    outfile = pysam.AlignmentFile(outfile_path, 'wb', template=samfile)
    final_reads = []

    itr = samfile.fetch(refname, int(min_site), int(max_site)+1)

    for x in itr:
        if x.query_name in reads_considered:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)

    samfile.close()
    outfile.close()

    print(f'Output saved to {outfile_path}')

    return final_reads
