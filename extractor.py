import os

from Bio import GenBank
from Bio.Seq import Seq

from Feature import Feature
from main import dir_legalize, dir_with_time
from statistictor import write_out_file
from Bio.SeqUtils import GC


def extract(filename: str, contain: list = None, exceptl: list = None):
    # contain：只包含这些类型 contain annotations in this list only
    # exceptl：不包含这些类型 do not contain annotations in this list
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    root = os.path.join('output', dir_with_time(os.path.split(os.path.splitext(filename)[0])[1]) )
    root = dir_legalize(root)
    # 输出文件夹：/output/输入文件名+时间
    # output dir: /output/{input filename}-Time
    if not os.path.exists(root):
        os.makedirs(root)

    with open(filename) as file:
        seq_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为注释名称，value为出现次数
        # keys are sequence name, values are child dictionaries
        # keys of child dictionaries are annotation name, values are occurrence times
        anno_list: list[str] = []  # 注释名列表，内容为注释名称 items are all annotations
        info_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为序列属性，value为属性值
        # keys are sequence name, values are child dictionaries, which contains other properties and values
        info_list: list[str] = ['Length', 'LSC', 'SSC', 'IR', 'GC-content', 'PCG', 'tRNA', 'rRNA']

        for record in GenBank.parse(file):
            processed = []
            # 储存已处理过的feature，用于判断拷贝数，每条序列开始时重置。
            # processed feature, used for copy amount, reset when new sequence starts

            seq_name = (record.organism + '-' + record.accession[0]).replace(' ', '_')
            print(f'Parsing: {seq_name}')

            if seq_name not in seq_dict:
                seq_dict[seq_name] = {}
            if seq_name not in info_dict:
                info_dict[seq_name] = {'Length': '', 'LSC': '', 'SSC': '', 'IR': '', 'GC-content': '',
                                       'PCG': 0, 'tRNA': 0, 'rRNA': 0}
            info_dict[seq_name]['Length'] = len(record.sequence)
            info_dict[seq_name]['GC-content'] = GC(record.sequence)

            for feature in record.features:
                # 遇到list为空也会算作not in contain，所以应该先检查contain非空。
                if (len(contain) > 0 and feature.key.casefold() not in (ctype.casefold() for ctype in contain)) or (
                        feature.key.casefold() in (etype.casefold() for etype in exceptl)):
                    continue

                try:
                    myfeature = Feature(feature)
                    # 名称由基因名+注释类型组成，例如 ycf2 gene, ycf2 CDS
                    # fullname is gene name + annotation type, i. e. ycf2 gene, ycf2 CDS
                    if 'gene' in myfeature.qualifier_dict:
                        name_key = 'gene'
                    elif 'product' in myfeature.qualifier_dict:
                        name_key = 'product'
                    elif 'organism' in myfeature.qualifier_dict:
                        name_key = 'organism'
                    elif 'note' in myfeature.qualifier_dict:
                        name_key = 'note'
                    elif 'rpt_type' in myfeature.qualifier_dict:
                        name_key = 'rpt_type'

                    fullname = myfeature.qualifier_dict[name_key] + ' ' + myfeature.key

                    # print(f'Parsing: {fullname}')

                    if fullname not in anno_list:
                        anno_list.append(fullname)

                    # 该注释在该序列中出现的次数
                    # occurrence times of this annotation in this sequence
                    if fullname not in seq_dict[seq_name]:
                        seq_dict[seq_name][fullname] = 1
                    else:
                        seq_dict[seq_name][fullname] += 1

                    suffix = ''  # 后缀默认为空字符串
                    # print(record.organism + ' ' + record.accession[0] + ' ' + fullname + ' ', end=' ')

                    i = 2
                    while (fullname + suffix).casefold() in (item.casefold() for item in processed):
                        suffix = '_copy' + str(i)
                        i = i + 1

                    # for interval in myfeature.intervals:
                    # print(interval, end=' ')
                    # print()

                    # /output/filename-time/feature-type/gene-name.fasta
                    #   |--CDS
                    #     |--accD
                    #     |--rps12
                    #     |--ycf1
                    #     |--ycf2
                    #     |--ycf2_copy2
                    #     |--others
                    #   |--gene
                    #   |--rRNA
                    #   |--tRNA
                    #   |--others

                    keydir = dir_legalize(os.path.join(root, myfeature.key))
                    if not os.path.exists(keydir):
                        os.makedirs(keydir)
                    outname = os.path.join(keydir, myfeature.qualifier_dict[name_key] + suffix + '.fasta')
                    outname = dir_legalize(outname)

                    with open(outname, 'a', newline='') as outfile:  # 'w': write,写入并覆盖原有内容；‘a': append, 附加在原有内容之后。
                        outfile.write(f'>{seq_name}\n')
                        seq_of_it = ''
                        for interval in myfeature.intervals:
                            if interval.complement:
                                # 计算已知序列的反向互补序列 reverse complement of the sequence
                                seq_of_it = seq_of_it + str(
                                    Seq(record.sequence[interval.start - 1:interval.end]).reverse_complement())
                            else:
                                seq_of_it = seq_of_it + record.sequence[interval.start - 1:interval.end]

                        if 'LSC' in outname:
                            info_dict[seq_name]['LSC'] = len(seq_of_it)
                        elif 'SSC' in outname:
                            info_dict[seq_name]['SSC'] = len(seq_of_it)
                        elif 'IR' in outname or 'repeat_region' in outname:
                            if info_dict[seq_name]['IR'] == '':
                                info_dict[seq_name]['IR'] = len(seq_of_it)
                            elif info_dict[seq_name]['IR'] != len(seq_of_it):
                                print('repeat region length different, please check it')
                        elif 'CDS' in outname:
                            info_dict[seq_name]['PCG'] = info_dict[seq_name]['PCG'] + 1
                        elif 'tRNA' in outname:
                            info_dict[seq_name]['tRNA'] = info_dict[seq_name]['tRNA'] + 1
                        elif 'rRNA' in outname:
                            info_dict[seq_name]['rRNA'] = info_dict[seq_name]['rRNA'] + 1

                        outfile.write(f'{seq_of_it}\n')

                    processed.append(fullname + suffix)

                except KeyError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))
                except IndexError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))

        print('\nWriting statistic csv file')
        s_outname: str = os.path.join(root, 'annotation_statistic.csv')
        write_out_file(s_outname, anno_list, seq_dict)

        print('\nChecking properties and writing to csv file')
        p_outname: str = os.path.join(root, 'properties_statistic.csv')
        write_out_file(p_outname, info_list, info_dict)

        print('\nAnnotation extraction complete.')


if __name__ == '__main__':
    extract("input.gb", contain=[
        # 'CDS', 'gene', 'tRNA', 'rRNA'
    ], exceptl=[
        'exon', 'intron'
    ])
