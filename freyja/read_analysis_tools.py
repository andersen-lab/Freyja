import re
import pysam
import gffpandas.gffpandas as gffpd
from Bio.Seq import Seq

def extract(query_mutations, input_bam, output, refname, same_read):
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
    print("Extracting read pairs with specified mutations")
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

            if x.cigarstring is None:
                # checks for a possible fail case
                continue
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
                            reads_considered.append(x)
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
                            reads_considered.append(x)
                            continue
                        i += int(m[0])

            c_dict = dict(zip(x.get_reference_positions(), seq))
            for s in sites_in:
                if snp_dict[s] == c_dict[s]:
                    reads_considered.append(x)
                    break
    samfile.close()

    # Same process, this time removing reads that don't contain all
    # query mutations.

    if same_read:
        temp = reads_considered
        reads_considered = []

        for x in temp:
            ref_pos = set(x.get_reference_positions())
            start = x.reference_start
            sites_in = list(ref_pos & set(snp_sites))
            seq = x.query_alignment_sequence
            cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)

            for ins in insertions:
                i = 0
                del_spacing = 0
                ins_spacing = 0
                insert_found = False
                for m in cigar:
                    if m[1] == 'I':
                        if (start+i+del_spacing-ins_spacing,
                                seq[i:i+int(m[0])]) == ins:
                            insert_found = True
                            break
                        i += int(m[0])
                        ins_spacing += int(m[0])
                    elif m[1] == 'D':
                        del_spacing += int(m[0])
                    elif m[1] == 'M':
                        i += int(m[0])
                if not insert_found:
                    break
            if len(insertions) > 0 and not insert_found:
                continue

            for _del in deletions:
                deletion_found = False
                i = x.reference_start
                for m in cigar:
                    if m[1] == 'M':
                        i += int(m[0])
                    elif m[1] == 'D':
                        if (i, int(m[0])) == _del:
                            deletion_found = True
                            break
                        i += int(m[0])
                if not deletion_found:
                    break
            if len(deletions) > 0 and not deletion_found:
                continue

            all_snps_found = True
            c_dict = dict(zip(x.get_reference_positions(), seq))
            for s in snp_sites:
                if s not in sites_in or snp_dict[s] != c_dict[s]:
                    all_snps_found = False
                    break
            if not all_snps_found:
                continue

            reads_considered.append(x)

    # Run again, this time also getting the paired read
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    outfile = pysam.AlignmentFile(output, 'wb', template=samfile)
    final_reads = []
    query_names = [x.query_name for x in reads_considered]

    itr = samfile.fetch(refname, min(all_sites), max(all_sites)+1)
    for x in itr:
        if x.query_name not in query_names:
            continue
        outfile.write(x)
        final_reads.append(x.query_name)

    samfile.close()
    outfile.close()

    print(f'Output saved to {output}')

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
        print('extract: Missing index file, try running samtools index',
              input_bam)
        return -1
    print("Filtering out reads with specified mutations")

    all_sites = snp_sites + indel_sites
    all_sites.sort()

    for site in all_sites:
        itr = samfile.fetch(refname, site, site+1)

        for x in itr:
            ref_pos = set(x.get_reference_positions())
            start = x.reference_start
            sites_in = list(ref_pos & set(snp_sites))

            seq = x.query_alignment_sequence

            if x.cigarstring is None:
                # checks for a possible fail case
                continue
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

    # performance-wise since we have to look at ALL reads.
    # Probably better to use a similar approach to extract
    # itr = samfile.fetch(refname, int(min_site), int(max_site)+1)

    # for x in itr:
    #     ref_pos = set(x.get_reference_positions())
    #     start = x.reference_start
    #     sites_in = list(ref_pos & set(snp_sites))
    #     if x.cigarstring is None:
    #         #checks for a possible fail case
    #         continue
    #     seq = x.query_alignment_sequence
    #     cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)

    #     # note: inserts add to read length, but not alignment coordinate
    #     if 'I' in x.cigarstring:
    #         i = 0
    #         del_spacing = 0
    #         ins_spacing = 0
    #         ranges = []
    #         insert_found = False
    #         for m in cigar:
    #             if m[1] == 'I':
    #                 if (start+i+del_spacing-ins_spacing,
    #                         seq[i:i+int(m[0])]) in insertions:
    #                     reads_considered.append(x.query_name)
    #                     insert_found = True
    #                     break
    #                 i += int(m[0])
    #                 ins_spacing += int(m[0])
    #             elif m[1] == 'D':
    #                 del_spacing += int(m[0])
    #                 continue
    #             elif m[1] == 'M':
    #                 ranges.append([i, i+int(m[0])])
    #                 i += int(m[0])

    #         seq = ''.join([seq[r[0]:r[1]] for r in ranges])

    #         if insert_found:
    #             continue

    #     if 'D' in x.cigarstring:
    #         i = x.reference_start
    #         for m in cigar:
    #             if m[1] == 'M':
    #                 i += int(m[0])
    #             elif m[1] == 'D':
    #                 if (i, int(m[0])) in deletions:
    #                     reads_considered.append(x.query_name)
    #                     continue
    #                 i += int(m[0])
    #     else:
    #         c_dict = dict(zip(x.get_reference_positions(), seq))
    #         for s in sites_in:
    #             if snp_dict[s] == c_dict[s]:
    #                 reads_considered.append(x.query_name)
    #                 break
    # samfile.close()

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

    print(f'Output saved to {output}')

    return final_reads


def get_cooccurrences(input_bam, min_site, max_site, output, refname, gff_file):
    
    if gff_file is not None:
        gene_positions = {}
        gff = gffpd.read_gff3(gff_file).filter_feature_of_type(['gene'])
        for i, attr in enumerate(gff.df['attributes']):
            attr = attr.split(';')
            gene_name = attr[2][5:]
            start_pos = gff.df['start'].iloc[i]
            end_pos = gff.df['end'].iloc[i]

            gene_positions[gene_name] = (start_pos, end_pos)
        
        if refname == 'NC_045512.2' and 'ORF1ab' in gene_positions:
            del gene_positions['ORF1ab']
            gene_positions['ORF1a'] = (266, 13468) 
            gene_positions['ORF1b'] = (13468, 21555)

    print(gene_positions)
    def get_gene(locus):
        for gene in gene_positions:
            start, end = gene_positions[gene]
            if locus in range(start, end+1):
                return gene, start

    try:
        samfile = pysam.AlignmentFile(input_bam, 'rb')
    except ValueError:
        print('get_cooccurrences: Missing index file. Try running samtools\
              index',
              input_bam)
        return -1
    


    co_muts = {}
    itr = samfile.fetch(refname, min_site, max_site+1)
    for x in itr:
        muts_found = []
        ref_seq = x.get_reference_sequence()
        start = x.reference_start
        seq = x.query_alignment_sequence

        if x.cigarstring is None:
            # checks for a possible fail case
            continue
        cigar = re.findall(r'(\d+)([A-Z]{1})', x.cigarstring)

        if 'I' in x.cigarstring:
                i = 0
                del_spacing = 0
                ins_spacing = 0
                ranges = []
    
                for m in cigar:
                    if m[1] == 'I':
                        muts_found.append(
                            (start+i+del_spacing-ins_spacing,
                             seq[i:i+int(m[0])])
                        )
                        i += int(m[0])
                        ins_spacing += int(m[0])
                    elif m[1] == 'D':
                        del_spacing += int(m[0])
                        continue
                    elif m[1] == 'M':
                        ranges.append([i, i+int(m[0])])
                        i += int(m[0])

                seq = ''.join([seq[r[0]:r[1]] for r in ranges])
        if 'D' in x.cigarstring:
                i = 0
                for m in cigar:
                    if m[1] == 'M':
                        i += int(m[0])
                    elif m[1] == 'D':
                        muts_found.append((start+i, int(m[0])))
                        i += int(m[0])
        
        seq = x.query_alignment_sequence
        # Find SNPs
        if ('S' not in x.cigarstring): # Ignore softclipped reads for now
            for tup in x.get_aligned_pairs(matches_only=True, with_seq=True):
                
                read_site, ref_site, ref_base = tup
                if ref_base != seq[read_site]:
                    muts_found.append(f'{ref_base.upper()}{ref_site+1}{seq[read_site]}')

        # Add gene-aa info to muts
        ref_seq  = Seq(x.get_reference_sequence())
        if gff_file is not None:
            for mut in muts_found:
                if type(mut) == tuple: 
                    if any(n in str(mut) for n in 'ACGT'): # Insertion
                        pass
                    else: # Deletion
                        pass
                else: # SNP
                    locus = int(mut[1:-1])
                    gene_info = get_gene(locus)
                    if gene_info is None: # locus occurs outside a gene
                        continue

                    gene, start_site = gene_info


        name = ','.join([str(mut) for mut in muts_found])
        if len(name) > 1:
            if name not in co_muts:
                co_muts[name] = 1
            else:
                co_muts[name] += 1
            

    thresh = 250
    for k in co_muts:
        if co_muts[k] > thresh:
            print(k, co_muts[k])
    samfile.close()