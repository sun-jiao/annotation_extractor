import os
import time

from Bio import SeqIO


def combine(infile_list: list, molecular_type: str):
    seq_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为基因名称，value为碱基序列
    # keys are sequence name, values are child dictionaries
    # keys of child dictionaries are gene(or CDS, trnA, etc) name, values are sequence
    gene_dict: dict = {}  # key为基因名称，value为该基因的长度
    # keys are all listed genes, values are longest length of the gene

    for infile in infile_list:
        for fasta_item in os.listdir(infile):
            if not (fasta_item.endswith('.fasta') or fasta_item.endswith('.fas') or fasta_item.endswith('.fsa')):
                continue

            gene_name = fasta_item.split('.')[0]
            if gene_name not in gene_dict:
                gene_dict[gene_name] = 0
            for seq_record in SeqIO.parse(os.path.join(infile, fasta_item), "fasta"):
                if seq_record.id not in seq_dict:
                    seq_dict[seq_record.id] = {}
                seq_dict[seq_record.id][gene_name] = str(seq_record.seq)
                if len(seq_record) > gene_dict[gene_name]:
                    gene_dict[gene_name] = len(seq_record)

    all_dict: dict = {}  # key 为序列名称，value为合并后的完整的序列
    name_length = 0
    for seq_name in seq_dict:
        child_dict = seq_dict[seq_name]
        if len(seq_name) > name_length:
            name_length = len(seq_name)
        if seq_name not in all_dict:
            all_dict[seq_name] = ''
        for gene_name in gene_dict:
            length = gene_dict[gene_name]
            if gene_name in child_dict:
                seq = child_dict[gene_name]
            else:
                seq = ''

            while len(seq) < length:
                seq = seq + '-'
            all_dict[seq_name] = all_dict[seq_name] + seq

    file_dir = os.path.join('output', 'combine-' + str(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())))
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)

    total_length = 0
    # writing partition file:
    part_filename = os.path.join(file_dir, 'partition.csv')

    with open(part_filename, 'w', newline='') as outfile:
        outfile.write('gene name,length,start,end\n')
        for gene_name in gene_dict:
            length = gene_dict[gene_name]
            outfile.write(f'{gene_name},{length},{total_length + 1},')
            total_length = total_length + length
            outfile.write(f'{total_length}\n')

    # writing fasta file:
    fasta_filename = os.path.join(file_dir, 'combine.fasta')

    with open(fasta_filename, 'w', newline='') as outfile:
        for name in all_dict:
            outfile.write(f'>{name_format(name)}\n')
            outfile.write(f'{all_dict[name]}\n')

    # writing nex file:
    nex_filename = os.path.join(file_dir, 'combine.nex')

    with open(nex_filename, 'w', newline='') as outfile:
        outfile.write('#' + f'NEXUS\nBEGIN DATA;\ndimensions ntax={str(len(all_dict))} nchar={str(total_length)};\n'
                            f'format missing=?\ndatatype={molecular_type} gap= -;\n\nmatrix\n')
        for name in all_dict:
            key = name
            name = name_format(name)
            while len(name) < name_length + 1:
                name = name + ' '
            outfile.write(f'{name}{all_dict[key]}\n')
        outfile.write(';\nEND;\n\n')

    # writing phy file:
    phy_filename = os.path.join(file_dir, 'combine.phy')

    with open(phy_filename, 'w', newline='') as outfile:
        outfile.write(f' {str(len(all_dict))} {str(total_length)}\n')
        for name in all_dict:
            key = name
            name = name_format(name)
            while len(name) < name_length + 1:
                name = name + ' '
            outfile.write(f'{name}{all_dict[key]}\n')
        outfile.write('\n')


def name_format(name: str):
    name = name.replace('Unclassified.', '_').replace('[^A-Za-z_]', '_')
    while '__' in name:
        name = name.replace('__', '_')
    return name


if __name__ == '__main__':
    combine(['./output/input-2021-07-15_13-55-37/CDS_MUSCLE'], 'DNA')
