import os
import time

from Bio import GenBank
from Bio.Seq import Seq

from Feature import Feature
from statistictor import write_out_file


def extract(filename: str, contain: list = None, exceptl: list = None):
    # contain：只包含这些类型 contain annotations in this list only
    # exceptl：不包含这些类型 do not contain annotations in this list
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    root = 'output/' + filename.split('.')[-2] + '-' + str(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()))
    # 输出文件夹：/output/输入文件名+时间
    # output dir: /output/{input filename}-Time
    if not os.path.exists(root):
        os.mkdir(root)

    with open(filename) as file:
        seq_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为注释名称，value为出现次数
        # keys are sequence name, values are child dictionaries
        # keys of hild dictionaries are annotation name, values are occurrence times
        anno_list: list[str] = []  # 注释名列表，内容为注释名称 items are all annotations

        for record in GenBank.parse(file):
            processed = []
            # 储存已处理过的feature，用于判断拷贝数，每条序列开始时重置。
            # processed feature, used for copy amount, reset when new sequence starts

            seq = (record.organism + '-' + record.accession[0]).replace(' ', '_')
            if seq not in seq_dict:
                seq_dict[seq] = {}

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

                    fullname = myfeature.qualifier_dict[name_key] + ' ' + myfeature.key

                    print(f'Parsing: {fullname}')

                    if fullname not in anno_list:
                        anno_list.append(fullname)

                    # 该注释在该序列中出现的次数
                    # occurrence times of this annotation in this sequence
                    if fullname not in seq_dict[seq]:
                        seq_dict[seq][fullname] = 1
                    else:
                        seq_dict[seq][fullname] += 1

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

                    outname = root + '\\' + myfeature.key + '\\' + myfeature.qualifier_dict[name_key] + suffix + '.fasta'
                    if not os.path.exists(root + '\\' + myfeature.key):
                        os.mkdir(root + '\\' + myfeature.key)

                    with open(outname, 'a', newline='') as outfile:  # 'w': write,写入并覆盖原有内容；‘a': append, 附加在原有内容之后。
                        outfile.write('>' + seq + '\n')
                        for interval in myfeature.intervals:
                            if interval.complement:
                                outfile.write(
                                    # 计算已知序列的反向互补序列 reverse complement of the sequence
                                    str(Seq(record.sequence[interval.start - 1:interval.end]).reverse_complement()))
                            else:
                                outfile.write(record.sequence[interval.start - 1:interval.end])
                        outfile.write('\n')

                    processed.append(fullname + suffix)

                except KeyError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))
                except IndexError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))

        print('Writing csv file')
        outname: str = root + '/annotation_statistic.csv'

        write_out_file(outname, anno_list, seq_dict)
        print('Annotation extraction complete.')


if __name__ == '__main__':
    extract("input.gb", contain=[# 'CDS', 'gene', 'tRNA', 'rRNA'
                                 ], exceptl=[])
