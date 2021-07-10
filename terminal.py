#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys, getopt

import extractor
import macse_align
import statistictor


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hp:i:c:e:t:h:",
                                   ["program=", "infile=", "contain=", "except=", "transl_table=", "help="])
    except getopt.GetoptError:
        args_error()
        return

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('python terminal.py -p/--program program_name [args]\r\n'
                  'program list:\r\n'
                  '\textract\r\n'
                  '\tstatistic\r\n'
                  '\tmacse\r\n'
                  'extract and statistic args:\r\n'
                  '\t-i --infile\tinput file\r\n'
                  '\t-c --contain\tcontain annotations in this list only, separate using \',\', do not use space\r\n'
                  '\t-e --except\tdo not contain annotations in this list\r\n'
                  'macse args:\r\n'
                  '\t-i --infile\tinput file\r\n'
                  '\t-t --transl_table\ttranslation table (listed below)\r\n'
                  '\t\t1\tThe_Standard_Code\r\n'
                  '\t\t2\tThe_Vertebrate_Mitochondrial_Code\r\n'
                  '\t\t3\tThe_Yeast_Mitochondrial_Code\r\n'
                  '\t\t4\tThe_Mold_Protozoan_and_Coelenterate_Mitochondrial_Code_and_the_Mycoplasma_Spiroplasma_Code\r\n'
                  '\t\t5\tThe_Invertebrate_Mitochondrial_Code\r\n'
                  '\t\t6\tThe_Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code\r\n'
                  '\t\t9\tThe_Echinoderm_and_Flatworm_Mitochondrial_Code\r\n'
                  '\t\t10\tThe_Euplotid_Nuclear_Code\r\n'
                  '\t\t11\tThe_Bacterial_Archaeal_and_Plant_Plastid_Code\r\n'
                  '\t\t12\tThe_Alternative_Yeast_Nuclear_Code\r\n'
                  '\t\t13\tThe_Ascidian_Mitochondrial_Code\r\n'
                  '\t\t14\tThe_Alternative_Flatworm_Mitochondrial_Code\r\n'
                  '\t\t15\tBlepharisma_Nuclear_Code\r\n'
                  '\t\t16\tChlorophycean_Mitochondrial_Code\r\n'
                  '\t\t21\tTrematode_Mitochondrial_Code\r\n'
                  '\t\t22\tScenedesmus_obliquus_mitochondrial_Code\r\n'
                  '\t\t23\tThraustochytrium_Mitochondrial_Code\r\n')
            sys.exit()
        elif opt in ("-p", "--program"):
            if arg in ('extract', 'statistic'):
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
                    if arg == 'extract':
                        extractor.extract(infile, contain, exceptl)
                    elif arg == 'statistic':
                        statistictor.csv_statistics(infile, contain, exceptl)
            elif arg == 'macse':
                infile = ''
                transl_table = ''
                for opt0, arg0 in opts:
                    if opt0 in ('-i', '--infile'):
                        infile = arg0
                    elif opt0 in ('-t', '--transl_table'):
                        transl_table = arg0
                if infile != '':
                    macse_align.macse(infile, transl_table)
            else:
                args_error()


def args_error():
    print('arguments error, use \'-h\' or \'--help\' for more information.')
    sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])
