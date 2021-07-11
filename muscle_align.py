import os


def muscle(infile: str):
    if not os.path.isfile('muscle.exe'):
        print('Please download latest version of muscle from: '
              'https://www.drive5.com/muscle/downloads.htm \n'
              'place it in annotation_extractor dir, and rename it to \"muscle.exe\".')
        return

    CDS_path = '.\\' + infile + '\\CDS\\'
    NUC_out = '.\\' + infile + '\\CDS_MUSCLE\\'

    if not os.path.exists(NUC_out):
        os.mkdir(NUC_out)

    amount = len(os.listdir(CDS_path))
    for index, filename in enumerate(os.listdir(CDS_path)):
        print(f'{index} of {amount} complete' + '({:.2%})'.format(index/amount) + f'\nAligning {filename}')
        os.system('.\muscle -in ' + CDS_path + filename
                  + ' -out ' + NUC_out + filename)
    print(f'{amount} of {amount} complete (100%)')


if __name__ == '__main__':
    muscle('output\\input-2021-07-10_21-19-07')

