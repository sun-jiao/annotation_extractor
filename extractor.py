import os
import time

from Bio import GenBank
from Bio.Seq import Seq

from Feature import Feature


def extract(filename: str, contain: list = None, exceptl: list = None, statistic: bool = True):
    # contain：只包含这些类型；exceptl：不包含这些类型；两个同时存在的，前者优先。
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    root = 'output-' + filename.split('.')[0] + '-' + str(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()))
    # 输出文件夹：输入文件名+时间
    if not os.path.exists(root):
        os.mkdir(root)

    with open(filename) as file:
        for record in GenBank.parse(file):
            processed = []  # 储存已处理过的feature，用于判断拷贝数，每条序列开始时重置。
            for feature in record.features:
                # 遇到list为空也会算作not in contain，所以应该先检查contain非空。
                if (len(contain) > 0 and feature.key.casefold() not in (ctype.casefold() for ctype in contain)) or (
                        feature.key.casefold() in (etype.casefold() for etype in exceptl)):
                    continue

                try:
                    myfeature = Feature(feature)
                    fullname = myfeature.qualifier_dict[
                                   'gene'] + ' ' + myfeature.key  # 名称由基因名+注释类型组成，例如 ycf2 gene, ycf2 CDS
                    suffix = ''  # 后缀默认为空字符串
                    # print(record.organism + ' ' + record.accession[0] + ' ' + fullname + ' ', end=' ')

                    i = 2
                    while (fullname + suffix).casefold() in (item.casefold() for item in processed):
                        suffix = '_copy' + str(i)
                        i = i + 1

                    # for interval in myfeature.intervals:
                    # print(interval, end=' ')
                    # print()

                    # 文件目录结构：/filename-time/feature-type/gene-name.fasta
                    #  --root
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

                    outname = root + '\\' + myfeature.key + '\\' + myfeature.qualifier_dict['gene'] + suffix + '.fasta'
                    if not os.path.exists(root + '\\' + myfeature.key):
                        os.mkdir(root + '\\' + myfeature.key)

                    with open(outname, 'a', newline='') as outfile:  # 'w': write,写入并覆盖原有内容；‘a': append, 附加在原有内容之后。
                        outfile.write('>' + record.organism + '-' + record.accession[0] + '\r\n')
                        for interval in myfeature.intervals:
                            if interval.complement:
                                outfile.write(
                                    # 计算已知序列的反向互补序列
                                    str(Seq(record.sequence[interval.start - 1:interval.end]).reverse_complement()))
                            else:
                                outfile.write(record.sequence[interval.start - 1:interval.end])
                        outfile.write('\r\n')

                    processed.append(fullname + suffix)

                except KeyError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\r\n' + str(feature))
                except IndexError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\r\n' + str(feature))


if __name__ == '__main__':
    extract("input.gb", contain=['CDS', 'gene', 'tRNA', 'rRNA'], exceptl=[])
