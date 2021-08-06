import os

from main import app_name


def macse(input: str, transl_table: str):
    if not os.path.isfile('macse.jar'):
        print('Please download latest version of macse from: '
              'https://bioweb.supagro.inra.fr/macse/index.php?menu=releases \n'
              f'place it in {app_name} dir, and rename it to \"macse.jar\".')
        return

    CDS_path = os.path.join(input, 'CDS')
    AA_out = os.path.join(input, 'CDS_AA')
    NUC_out = os.path.join(input, 'CDS_NUC')
    if (not os.path.exists(CDS_path)) or len(os.listdir(CDS_path)) == 0:
        print('CDS folder not exist or empty, nothing to align.')
        return
    if not os.path.exists(AA_out):
        os.makedirs(AA_out)
    if not os.path.exists(NUC_out):
        os.makedirs(NUC_out)

    amount = len(os.listdir(CDS_path))
    for index, filename in enumerate(os.listdir(CDS_path)):
        print(f'{index} of {amount} complete' + '({:.2%})'.format(index/amount))
        os.system('java -jar macse.jar -prog alignSequences -seq ' + CDS_path + filename
                  + ' -out_NT ' + NUC_out + filename + ' -out_AA ' + AA_out + filename + ' -gc_def ' + transl_table)
    print(f'{amount} of {amount} complete (100%)')


if __name__ == '__main__':
    macse('\\output\\input-2021-07-08_21-39-56', str(11))

