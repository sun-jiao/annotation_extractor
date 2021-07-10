import os
import subprocess


def macse(input: str, transl_table: str):
    if not os.path.isfile('macse.jar'):
        print('Please download latest version of macse from: '
              'https://bioweb.supagro.inra.fr/macse/index.php?menu=releases \n'
              'place it in annotation_extractor dir, and rename it to \"macse.jar\".')
        return

    CDS_path = '.\\' + input + '\\CDS\\'
    AA_out = '.\\' + input + '\\CDS_AA\\'
    NUC_out = '.\\' + input + '\\CDS_NUC\\'
    if not os.path.exists(AA_out):
        os.mkdir(AA_out)
    if not os.path.exists(NUC_out):
        os.mkdir(NUC_out)

    amount = len(os.listdir(CDS_path))
    for index, filename in enumerate(os.listdir(CDS_path)):
        print(f'{index} of {amount} complete' + '({:.2%})'.format(index/amount))
        os.system('java -jar macse.jar -prog alignSequences -seq ' + CDS_path + filename
                  + ' -out_NT ' + NUC_out + filename + ' -out_AA ' + AA_out + filename + ' -gc_def ' + transl_table)
    print(f'{amount} of {amount} complete (100%)')


if __name__ == '__main__':
    macse('input-2021-07-08_21-39-56', str(11))

