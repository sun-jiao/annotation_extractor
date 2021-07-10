import csv
import os


def csv_statistics(filename: str, contain: list = None, exceptl: list = None):
    # contain：只包含这些类型；exceptl：不包含这些类型；两个同时存在的，前者优先。
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    with open(filename, newline='') as csvfile:
        seq_dict: dict = {}  # key为序列名称，value为另一个字典，子字典的key为注释名称，value为出现次数
        anno_list: list[str] = []  # 注释名列表，内容为注释名称

        csvreader = csv.reader(csvfile)
        for row in csvreader:
            # geneious导出的csv文件，每行第一个为Sequence Name，第二个为Name（Annotation Name），第三个为Type
            # 如果使用的是其它软件导出的，可自行修改此处的index
            seq: str = row[0]
            anno: str = row[1]
            type: str = row[2]

            if seq not in seq_dict:
                seq_dict[seq] = {}

            # 遇到list为空也会算作not in contain，所以应该先检查contain非空。
            if (len(contain) > 0 and type.casefold() not in (ctype.casefold() for ctype in contain)) or (
                    type.casefold() in (etype.casefold() for etype in exceptl)):
                continue

            if anno not in anno_list:
                anno_list.append(anno)

            # 该注释在该序列中出现的次数
            if anno not in seq_dict[seq]:
                seq_dict[seq][anno] = 1
            else:
                seq_dict[seq][anno] += 1

        print("seq length: ", str(len(seq_dict)), "; anno length: ", str(len(anno_list)))

    outname: str = filename.replace(".csv", "_out.csv")

    write_out_file(outname, anno_list, seq_dict)


def write_out_file(outname: str, anno_list: list, seq_dict: dict):
    if not os.path.exists(outname):
        with open(outname, 'w'): pass

    with open(outname, 'w', newline='') as csvfile:  # 'w'表示写入
        writer = csv.writer(csvfile)
        csvfile_header: list = ['seq_name']
        csvfile_header.extend(anno_list)
        writer.writerow(csvfile_header)
        for key in seq_dict:
            value = seq_dict[key]
            csvfile_line: list = [key]
            for anno_name in anno_list:
                if anno_name in value:
                    csvfile_line.append(value[anno_name])
                else:
                    csvfile_line.append(0)
            writer.writerow(csvfile_line)


if __name__ == '__main__':
    csv_statistics('input.csv', contain=[], exceptl=['intron', 'exon', 'misc_feature', 'repeat_region'])
