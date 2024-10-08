# Annotation Extractor

This is a simple biology tool for extract all annotations you want from Genbank file(*.gb). Based on Biopython package.

## Usage

basic usage:

    python terminal.py -p/--program program_name [args]

arguments:

    program list:
        extract: extract all annotations you want to fasta files
        statistic: statistic occurence times of annotations
        muscle: run muscle alignment
        macse: run macse alignment

    extract, statistic and muscle args:
        -i --infile     input file
        -c --contain    contain annotations in this list only, separate using ',', do not use space
        -e --except     do not contain annotations in this list

    combine args:
        -i --infile     input file folders list, separate using ',', do not use space
        -m --mole_type  molecular type, should be one of DNA, RNA, Protein

    macse args:
        -i --infile     input file
        -t --transl_table       translation table (listed below)
                1       The_Standard_Code
                2       The_Vertebrate_Mitochondrial_Code
                3       The_Yeast_Mitochondrial_Code
                4       The_Mold_Protozoan_and_Coelenterate_Mitochondrial_Code_and_the_Mycoplasma_Spiroplasma_Code
                5       The_Invertebrate_Mitochondrial_Code
                6       The_Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code
                9       The_Echinoderm_and_Flatworm_Mitochondrial_Code
                10      The_Euplotid_Nuclear_Code
                11      The_Bacterial_Archaeal_and_Plant_Plastid_Code
                12      The_Alternative_Yeast_Nuclear_Code
                13      The_Ascidian_Mitochondrial_Code
                14      The_Alternative_Flatworm_Mitochondrial_Code
                15      Blepharisma_Nuclear_Code
                16      Chlorophycean_Mitochondrial_Code
                21      Trematode_Mitochondrial_Code
                22      Scenedesmus_obliquus_mitochondrial_Code
                23      Thraustochytrium_Mitochondrial_Code

### Examples

extractor, extract annotations from input.gb, contain CDS, gene, rRNA, tRNA, except repeat_region and misc_feature:

    python terminal.py -p extract -i input.gb -c CDS,gene,rRNA,tRNA -e repeat_region,misc_feature

macse alignment, translate table is 11:

    python terminal.py -p macse -i \output\input-20xx-mm-dd_hh-mm-ss -t 11

combine all fasta files in CDS folder:

    python terminal.py -p combine -i \CDS,\tRNA -m DNA

## Others

Most of untranslated Chinese comments are not important. Please don't care about it.
