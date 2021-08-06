import os
import re
import time

app_name = 'annotation extractor'


# os.makedirs() 递归创建目录

def dir_legalize(filename: str):
    dir = re.split(r'[/\\]', filename)
    newdir = []
    for level in dir:
        newdir.append(legalize(level))
    return os.path.join(*newdir)


def legalize(string: str):
    return re.sub(r'[^A-Za-z0-9_]', '_', string)


def dir_with_time(prefix: str):
    return prefix + '_' + str(time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()))
