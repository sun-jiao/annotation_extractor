#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys, getopt

import extractor
import macse_align
import muscle_align
import statistictor
from conbine import combine
from main import app_name


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hp:i:c:e:t:h:m:",
                                   ["program=", "infile=", "contain=", "except=", "transl_table=", "help=", "mole_type="])
    except getopt.GetoptError:
        args_error()
        return

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(f'{app_name} usage help:\n'
                  'python terminal.py -p/--program program_name [args]\n'
                  'program list:\n'
                  '\textract: extract all annotations you want to fasta files\n'
                  '\tstatistic: statistic occurence times of annotations\n'
                  '\tmuscle: run muscle alignment\n'
                  '\tmacse: run macse alignment\n'
                  '\tcombine: combine multiple fasta file into fasta, nexus, phy files\n'
                  'extract, statistic and muscle args:\n'
                  '\t-i --infile\tinput file\n'
                  '\t-c --contain\tcontain annotations in this list only, separate using \',\', do not use space\n'
                  '\t-e --except\tdo not contain annotations in this list\n'
                  'combine args:\n'
                  '\t-i --infile\tinput file folders list, separate using \',\', do not use space\n'
                  '\t-m --mole_type\tmolecular type, should be one of DNA, RNA, Protein\n'
                  'macse args:\n'
                  '\t-i --infile\tinput file\n'
                  '\t-t --transl_table\ttranslation table (listed below)\n'
                  '\t\t1\tThe_Standard_Code\n'
                  '\t\t2\tThe_Vertebrate_Mitochondrial_Code\n'
                  '\t\t3\tThe_Yeast_Mitochondrial_Code\n'
                  '\t\t4\tThe_Mold_Protozoan_and_Coelenterate_Mitochondrial_Code_and_the_Mycoplasma_Spiroplasma_Code\n'
                  '\t\t5\tThe_Invertebrate_Mitochondrial_Code\n'
                  '\t\t6\tThe_Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code\n'
                  '\t\t9\tThe_Echinoderm_and_Flatworm_Mitochondrial_Code\n'
                  '\t\t10\tThe_Euplotid_Nuclear_Code\n'
                  '\t\t11\tThe_Bacterial_Archaeal_and_Plant_Plastid_Code\n'
                  '\t\t12\tThe_Alternative_Yeast_Nuclear_Code\n'
                  '\t\t13\tThe_Ascidian_Mitochondrial_Code\n'
                  '\t\t14\tThe_Alternative_Flatworm_Mitochondrial_Code\n'
                  '\t\t15\tBlepharisma_Nuclear_Code\n'
                  '\t\t16\tChlorophycean_Mitochondrial_Code\n'
                  '\t\t21\tTrematode_Mitochondrial_Code\n'
                  '\t\t22\tScenedesmus_obliquus_mitochondrial_Code\n'
                  '\t\t23\tThraustochytrium_Mitochondrial_Code\n')
            sys.exit()
        elif opt in ("-p", "--program"):
            if arg in ('extract', 'statistic', 'muscle'):
                infile = ''
                contain = ''
                exceptl = ''
                for opt0, arg0 in opts:
                    if opt0 in ('-i', '--infile'):
                        infile = arg0
                    elif opt0 in ('-c', '--contain'):
                        contain = arg0.split(',')
                    elif opt0 in ('-e', '--except'):
                        exceptl = arg0.split(',')
                if infile != '':
                    try:
                        if arg == 'extract':
                            extractor.extract(infile, contain, exceptl)
                        elif arg == 'statistic':
                            statistictor.csv_statistics(infile, contain, exceptl)
                        elif arg == 'muscle':
                            muscle_align.muscle(infile, contain, exceptl)
                    except ValueError:
                        args_error()
            elif arg == 'combine':
                infile = ''
                mole_type = ''
                for opt0, arg0 in opts:
                    if opt0 in ('-i', '--infile'):
                        infile = arg0.split(',')
                    elif opt0 in ('-m', '--mole_type'):
                        mole_type = str(arg0)
                if infile != '':
                    try:
                        combine(infile, mole_type)
                    except ValueError as e:
                        print(e)
                else:
                    args_error()
            elif arg == 'macse':
                infile = ''
                transl_table = ''
                for opt0, arg0 in opts:
                    if opt0 in ('-i', '--infile'):
                        infile = arg0
                    elif opt0 in ('-t', '--transl_table'):
                        transl_table = arg0
                if infile != '':
                    try:
                        macse_align.macse(infile, transl_table)
                    except ValueError:
                        args_error()
                else:
                    args_error()
            else:
                args_error()


def args_error():
    print('arguments error, use \'-h\' or \'--help\' for more information.')
    sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])
