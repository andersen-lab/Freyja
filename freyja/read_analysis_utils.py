import pandas as pd
from Bio.Seq import MutableSeq
from Bio import SeqIO


def nt_position(x):
    if ',' in x:
        return int(x.split(',')[0][1:])
    return int(x.split('(')[0][1:-1])


def get_colnames_and_sites(unique_mutations, nt_muts):
    if nt_muts:
        nt_muts = True
        colnames = unique_mutations
        sites = {mut: nt_position(mut) for mut in unique_mutations}
    else:
        colnames = [mut.split(')(')[1][:-1] if 'DEL' in mut
                    else mut.split('(')[1][:-1] for mut in unique_mutations]
        sites = {mut.split(')(')[1][:-1] if 'DEL' in mut
                 else mut.split('(')[1][:-1]: nt_position(mut)
                 for mut in unique_mutations}
    colnames = sorted(list(set(colnames)), key=lambda x: sites[x])
    return colnames, sites


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


def filter_covariants_output(cluster, nt_muts, min_mutations):
    cluster_final = []

    if nt_muts:
        for variant in cluster:
            if ',' in variant:
                # Insertion
                if any([nt in variant for nt in 'ACGT']):
                    continue
                # Frameshift
                if int(variant.split(',')[1].split(')')[0]) % 3 != 0:
                    continue
            cluster_final.append(variant)
    else:
        for variant in cluster:
            if '*' in variant or 'INS' in variant:
                continue  # Ignore stop codons and insertions for now
            elif 'DEL' in variant \
                    and int(variant.split(")(")[0].split(",")[1]) % 3 != 0:
                continue  # Ignore frameshift mutations
            else:
                cluster_final.append(variant)

    if len(cluster_final) < min_mutations:
        return pd.NA
    else:
        return list(dict.fromkeys(cluster_final))


def parse_gff(gff_file):
    gene_positions = {}
    with open(gff_file) as f:
        for line in f.readlines():
            line = line.split('\t')
            if 'gene' in line:
                attrs = line[-1].split(';')
                for attr in attrs:
                    if attr.startswith('Name='):
                        gene_name = attr.split('=')[1]
                        gene_positions[gene_name] = (int(line[3]),
                                                     int(line[4]))
    return gene_positions


def translate_snps(snps, ref, gene_positions):

    # Load reference genome
    ref = MutableSeq(next(SeqIO.parse(ref, 'fasta')).seq)

    output = {snp: None for snp in snps}
    for snp in snps:
        locus = int(snp[1:-1])

        # Find the key in gene_positions that corresponds to the gene
        # containing the SNP
        gene_info = None
        for gene in gene_positions:
            start, end = gene_positions[gene]
            if start <= locus <= end:
                gene_info = gene, start
                break
        if gene_info is None:
            continue

        codon_position = (locus - start) % 3
        aa_locus = ((locus - codon_position - start) // 3) + 1

        ref_codon = ref[locus - codon_position - 1:
                        locus - codon_position + 2]
        ref_aa = ref_codon.translate()

        alt_codon = MutableSeq(str(ref_codon))
        alt_codon[codon_position] = snp[-1]

        if len(alt_codon) % 3 != 0 or len(alt_codon) == 0 \
                or 'N' in alt_codon:
            continue
        alt_aa = alt_codon.translate()

        if ref_aa == alt_aa or '*' in alt_aa:
            continue # Synonymous/stop codon mutation

        aa_mut = f'{gene}:{ref_aa}{aa_locus}{alt_aa}'
        output[snp] = aa_mut

    return output
