import os
import re
import pysam

def extract(query_mutations, input_bam, output):
    # Load data
    with open(query_mutations) as infile:
        mutations = infile.read().splitlines()

        if len(mutations) < 3:
            for i in range(3 - len(mutations)):
                mutations.append('')
        snps = [s.strip() for s in mutations[0] .split(',')]
        insertions = [s.strip() for s in mutations[1].split(',')]
        deletions = [s.strip() for s in mutations[2].split(',')]
    try:
        # parse tuples from indels
        if len(insertions[0]) > 0:
            insertions = [(int(s.split(':')[0][1:]),
                        s.split(':')[1][1:-2].strip('\''))
                        for s in insertions]
        if len(deletions[0]) > 0:
            deletions = [(int(s.split(':')[0][1:]), int(s.split(':')[1][:-1]))
                        for s in deletions]

        # get loci for all mutations
        snp_sites = [int(m[1:len(m)-1])-1 for m in snps if m]
        indel_sites = [s[0] for s in insertions if s] +\
                    [s[0] for s in deletions if s]
    except:
        print('extract: Error parsing',query_mutations)
        print('extract: See README for formatting requirements.')
        return -1

    snp_dict = {int(mut[1:len(mut)-1])-1 : mut[-1] for mut in snps if mut}

    samfile = pysam.AlignmentFile(input_bam, 'rb')

    reads_considered = []

    # TODO: Only fetch reads that contain mutations of interest
    # e.g. loop over samfile.fetch("NC_045512.2", mySite,mySite+1)
    # for each mySite in sites

    if len(snp_sites) == 0:
        snp_sites = [0]
    if len(indel_sites) == 0:
        indel_sites = [0]
    min_site = min(min(snp_sites), min(indel_sites))
    max_site = max(max(snp_sites), max(indel_sites))
    itr = samfile.fetch("NC_045512.2", min_site, max_site+1)

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
                    ranges.append([i,i+int(m[0])])
                    i += int(m[0])

            seq = ''.join([seq[r[0]:r[1]] for r in ranges])

            if insert_found:
                continue
                
        if 'D' in x.cigarstring:
            c_dict = dict(zip(x.get_reference_positions(), seq))
            i = x.reference_start

            for m in cigar:
                if m[1] == 'M':
                    i += int(m[0])
                elif m[1] == 'D':
                    for k in range(i,i+int(m[0])):
                            c_dict[k] = '-'
                    
                    if (i, int(m[0])) in deletions:
                        reads_considered.append(x.query_name)
                        continue

                    i += int(m[0])
            
            for s in sites_in:
                if c_dict[s] == snp_dict[s]:
                    reads_considered.append(x.query_name)
                    break
                    
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
    outfile = pysam.AlignmentFile(outfile_path,'wb',template=samfile)
    final_reads = []
    
    itr = samfile.fetch('NC_045512.2', min_site, max_site+1)

    for x in itr:
        if x.query_name not in reads_considered:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)
    
    samfile.close()
    outfile.close()

    print(f'Output saved to {outfile_path}')

    return final_reads

def filter(query_mutations, input_bam, min_site, max_site, output):

    # Load data
    with open(query_mutations) as infile:
        mutations = infile.read().splitlines()

        if len(mutations) < 3:
            for i in range(3 - len(mutations)):
                mutations.append('')
        snps = [s.strip() for s in mutations[0] .split(',')]
        insertions = [s.strip() for s in mutations[1].split(',')]
        deletions = [s.strip() for s in mutations[2].split(',')]
    try:
        # parse tuples from indels
        if len(insertions[0]) > 0:
            insertions = [(int(s.split(':')[0][1:]),
                        s.split(':')[1][1:-2].strip('\''))
                        for s in insertions]
        if len(deletions[0]) > 0:
            deletions = [(int(s.split(':')[0][1:]), int(s.split(':')[1][:-1]))
                        for s in deletions]

        # get loci for all mutations
        snp_sites = [int(m[1:len(m)-1])-1 for m in snps if m]
        indel_sites = [s[0] for s in insertions if s] +\
                    [s[0] for s in deletions if s]
    except:
        print('filter: Error parsing',query_mutations)
        print('filter: See README for formatting requirements.')
        return -1

    snp_dict = {int(mut[1:len(mut)-1])-1 : mut[-1] for mut in snps if mut}
    reads_considered = []

    samfile = pysam.AlignmentFile(input_bam, 'rb')
    itr = samfile.fetch("NC_045512.2", int(min_site), int(max_site)+1)

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
                    ranges.append([i,i+int(m[0])])
                    i += int(m[0])

            seq = ''.join([seq[r[0]:r[1]] for r in ranges])

            if insert_found:
                continue
                
        if 'D' in x.cigarstring:
            c_dict = dict(zip(x.get_reference_positions(), seq))
            i = x.reference_start

            for m in cigar:
                if m[1] == 'M':
                    i += int(m[0])
                elif m[1] == 'D':
                    for k in range(i,i+int(m[0])):
                            c_dict[k] = '-'
                    
                    if (i, int(m[0])) in deletions:
                        reads_considered.append(x.query_name)
                        continue

                    i += int(m[0])
            
            for s in sites_in:
                if c_dict[s] == snp_dict[s]:
                    reads_considered.append(x.query_name)
                    break
                    
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
    outfile = pysam.AlignmentFile(outfile_path,'wb',template=samfile)
    final_reads = []
    
    itr = samfile.fetch('NC_045512.2', int(min_site), int(max_site)+1)

    for x in itr:
        if x.query_name in reads_considered:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)
    
    samfile.close()
    outfile.close()

    print(f'Output saved to {outfile_path}')

    return final_reads