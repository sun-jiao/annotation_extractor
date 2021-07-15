import os

from Bio import SeqIO


def combine(infile: str, contain: list = None, exceptl: list = None):
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    seq_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为基因名称，value为碱基序列
    # keys are sequence name, values are child dictionaries
    # keys of child dictionaries are gene(or CDS, trnA, etc) name, values are sequence
    gene_dict: dict = {}  # key为基因名称，value为 items are all listed genes

    for diritem in os.listdir(infile):
        if (len(contain) > 0 and diritem.casefold() not in (ctype.casefold() for ctype in contain)) or (
                diritem.casefold() in (etype.casefold() for etype in exceptl)):
            continue
        for fastaitem in os.listdir(os.path.join(infile, diritem)):
            if not (fastaitem.endswith('.fasta') or fastaitem.endswith('.fas') or fastaitem.endswith('.fsa')):
                continue

            gene_name = fastaitem.split('.')[0]
            if gene_name not in gene_dict:
                gene_dict[gene_name] = 0
            for seq_record in SeqIO.parse(os.path.join(infile, diritem, fastaitem), "fasta"):
                if seq_record.id not in seq_dict:
                    seq_dict[seq_record.id] = {}
                seq_dict[seq_record.id][gene_name] = str(seq_record.seq)
                if len(seq_record) > gene_dict[gene_name]:
                    gene_dict[gene_name] = len(seq_record)

    all_dict: dict = {}
    for gene_name, length in gene_dict:
        for seq_name, child_dict in seq_dict:
            pass


if __name__ == '__main__':
    combine("./output/input-2021-07-15_13-55-37/", ['CDS_MUSCLE'], [])
