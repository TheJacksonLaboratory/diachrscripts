import argparse
import os
from collections import defaultdict
from csv import DictReader
from statistics import mean



class GeneToBait:
    def __init__(self, ensembl_gene_id, bait_id, values_dict) -> None:
        self._ensembl_gene_id = ensembl_gene_id
        self._bait_id = bait_id
        self._mean_expression_d = defaultdict(int)
        for celltype, expression_list in values_dict.items():
            self._mean_expression_d[celltype] = mean(expression_list)
            
    @property
    def ensembl_gene_id(self):
        return self._ensembl_gene_id
    
    @property
    def bait_id(self):
        return self._bait_id
    
    def get_expression_for_cell_type(self, celltype):
        if celltype not in self._mean_expression_d:
            raise ValueError(f"{celltype} not found")
        return self._mean_expression_d.get(celltype)
    
    def get_expression_string(self, celltypelist):
        """
        return a tab separated string with the values for the cell types
        """
        values = []
        for celltype in celltypelist:
            val = str(self.get_expression_for_cell_type(celltype))
            values.append(val)
        return "\t".join(values)



def parse_rnaseq(rnafile):
    n_bait_na = 0
    gene2bait_list = []
    tissues = {
        'Mac2.1': 'MAC2', 'Mac2.2': 'MAC2', 'Mac2.3': 'MAC2', 'Mac2.4': 'MAC2', 'Mac2.5': 'MAC2', 
        'Mon.1': 'MON', 'Mon.2': 'MON', 'Mon.3': 'MON', 'Mon.4': 'MON', 'Mon.5': 'MON', 
        'MK.1': 'MK', 'MK.2': 'MK', 'MK.3': 'MK', 
        'nCD4.1': 'CD4', 'nCD4.2': 'CD4', 'nCD4.3': 'CD4', 'nCD4.4': 'CD4', 'nCD4.5': 'CD4', 
        'nCD4.6': 'CD4', 'nCD4.7': 'CD4', 
        'nCD4.8': 'CD4', 'Ery.1': 'ERY', 'Ery.2': 'ERY', 
        'Mac1.1': 'MAC1', 'Mac1.2': 'MAC1', 'Mac1.3': 'MAC1', 'Mac1.4': 'MAC1', 
        'Mac0.1': 'MAC0', 'Mac0.2': 'MAC0', 'Mac0.3': 'MAC0', 'Mac0.4': 'MAC0', 
        'Neu.1': 'NEU', 'Neu.2': 'NEU', 'Neu.3': 'NEU', 'Neu.4': 'NEU', 'Neu.5': 'NEU', 
        'Neu.6': 'NEU',  'Neu.7': 'NEU'
    }
    with open(rnafile) as f:
        reader = DictReader(f, delimiter='\t')
        for row in reader:
            ensembl_gene_id = row['ENSEMBL_GENEID']
            bait_id = row['BaitID']
            if bait_id == 'NA':
                n_bait_na += 1
                continue
            values_d = defaultdict(list)
            for header, celltype in tissues.items():
                v = float(row[header])
                values_d[celltype].append(v)
            g2b = GeneToBait(ensembl_gene_id=ensembl_gene_id, bait_id=bait_id, values_dict=values_d)
            gene2bait_list.append(g2b)
    print(f"We extracted values for {len(gene2bait_list)} genes")
    return gene2bait_list


def parse_rmap(rmapfile):
    # This file has no header
    # chrom-start-end-digestID
    # define map with key: digestID and value: chr\tstart\tend
    digest_map = defaultdict(str)
    with open (rmapfile) as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(fields) != 4:
                raise ValueError(f"Malformed rmap line {line}")
            chrstring = "\t".join(fields[0:2])
            digestid = fields[3]
            digest_map[digestid] = chrstring
    print(f"We found {len(digest_map)} digests")
    return digest_map


def output_rna2bait_file(fname, gene2bait_list, digest_map):
    fh = open(fname, 'wt')
    celltypelist = ['MAC0', 'MAC1', 'MAC2', 'MON', 'MK', 'CD4', 'NEU']
    header = ["ensembl", "bait.id", "chr", "begin", "end"]
    header.extend(celltypelist)
    fh.write("\t".join(header) + "\n")
    for g2b in gene2bait_list:
        if g2b.bait_id not in digest_map:
            raise ValueError(f"Could not find bait {g2b.bait_id} in digest_map")
        chrstring = digest_map.get(g2b.bait_id)
        expression_string = g2b.get_expression_string(celltypelist)
        items = [g2b.ensembl_gene_id,  g2b.bait_id, chrstring, expression_string]
        fh.write("\t".join(items) + "\n")
    fh.close()
   

def main(rnafile, rmapfile):
    if not os.path.isfile(rnafile):
        raise ValueError("Must pass valid file as --rna argument")
    if not rnafile.endswith("GeneExpressionMatrix.txt"):
        raise ValueError("Must pass valid file as --rna /path/../GeneExpressionMatrix.txt")
    if not os.path.isfile(rmapfile):
        raise ValueError("Must pass valid file as --rmap argument")
    if not rmapfile.endswith("Digest_Human_HindIII.rmap"):
        raise ValueError("Must pass valid file as --rmap /path/../Digest_Human_HindIII.rmap")
    gene2bait_list = parse_rnaseq(rnafile) 
    digest_map = parse_rmap(rmapfile=rmapfile)
    fname = "rnaseq_to_bait.tsv"
    output_rna2bait_file(fname=fname, gene2bait_list=gene2bait_list, digest_map=digest_map)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rna', type=str, default="/home/peter/data/diachromatic/GeneExpressionMatrix.txt", required=False)
    parser.add_argument('--rmap', type=str, default="/home/peter/data/diachromatic/Human_hg19/Digest_Human_HindIII.rmap",required=False)
    args = parser.parse_args()
    main(args.rna, args.rmap)