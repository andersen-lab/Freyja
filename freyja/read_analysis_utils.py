import pandas as pd


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
