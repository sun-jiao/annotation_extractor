import os


def muscle(infile: str, contain: list = None, exceptl: list = None):
    if contain is None:
        contain = []
    if exceptl is None:
        exceptl = []

    if not os.path.isfile('muscle.exe'):
        print('Please download latest version of muscle from: '
              'https://www.drive5.com/muscle/downloads.htm \n'
              'place it in annotation_extractor dir, and rename it to \"muscle.exe\".')
        return

    for folder in os.listdir(infile):
        if os.path.isfile(os.path.join(infile , folder)) or folder.endswith('MUSCLE') or \
                (len(contain) > 0 and folder.casefold() not in (ctype.casefold() for ctype in contain)) or \
                (folder.casefold() in (etype.casefold() for etype in exceptl)):
            continue

        input_path = infile + '\\' + folder + '\\'
        out_path = '.\\' + infile + '\\' + folder + '_MUSCLE\\'

        if not os.path.exists(out_path):
            os.mkdir(out_path)

        print(f'Aligning {folder}:')

        amount = len(os.listdir(input_path))

        for index, filename in enumerate(os.listdir(input_path)):
            print(f'{index} of {amount} complete' + '({:.2%})'.format(index / amount) + f'\nAligning {filename}')
            os.system('.\muscle -in ' + input_path + filename
                      + ' -out ' + out_path + filename)
        print(f'{amount} of {amount} complete (100%)')


if __name__ == '__main__':
    muscle('output\\input-2021-07-10_21-19-07', contain=['CDS'])
