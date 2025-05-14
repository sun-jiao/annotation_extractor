# This script is for extracting the trnL-trnF marker.

import os

from Bio import GenBank
from typing import List

from Feature import Feature
from main import dir_legalize, dir_with_time
try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC


trnl = "trnL-UAA"
trnf = "trnF-GAA"

def extract(filename: str):
    root = os.path.join('output', dir_with_time(os.path.split(os.path.splitext(filename)[0])[1]) )
    root = dir_legalize(root)
    # 输出文件夹：/output/输入文件名+时间
    # output dir: /output/{input filename}-Time
    if not os.path.exists(root):
        os.makedirs(root)

    with open(filename) as file:
        for record in GenBank.parse(file):
            trnl_start = None
            trnf_start = None

            seq_name = (record.organism + '-' + record.accession[0]).replace(' ', '_')
            print(f'Parsing: {seq_name}')

            for feature in record.features:
                # 遇到list为空也会算作not in contain，所以应该先检查contain非空。

                try:
                    myfeature = Feature(feature)
                    # 名称由基因名+注释类型组成，例如 ycf2 gene, ycf2 CDS
                    # fullname is gene name + annotation type, i. e. ycf2 gene, ycf2 CDS
                    if myfeature.key.lower() != 'trna':
                        continue

                    if myfeature.qualifier_dict['gene'] == trnl:
                        trnl_start = myfeature.intervals[0].start
                    elif myfeature.qualifier_dict['gene'] == trnf:
                        trnf_start = myfeature.intervals[0].start
                    else:
                        continue

                    if trnl_start is None or trnf_start is None:
                        continue

                    outname = os.path.join(root, 'trnL_trnF.fasta')
                    outname = dir_legalize(outname)

                    with open(outname, 'a', newline='') as outfile:  # 'w': write,写入并覆盖原有内容；‘a': append, 附加在原有内容之后。
                        outfile.write(f'>{seq_name}\n')
                        seq_of_it = record.sequence[trnl_start - 1:trnf_start]
                        outfile.write(f'{seq_of_it}\n')

                    break

                except KeyError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))
                except IndexError as e:
                    print(str(e) + record.organism + ' ' + record.accession[0] + '\n' + str(feature))

        print('\nAnnotation extraction complete.')


if __name__ == '__main__':
    extract("127 documents.gb")
