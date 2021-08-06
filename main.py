import os
import re
import time

app_name = 'annotation extractor'


# os.makedirs() 递归创建目录
# os.path.split() 将路径和文件名二者分开，并非递归分割
# os.path.splitext() 将路径和扩展名二者分开

def dir_legalize(filename: str):
    splitext = os.path.splitext(filename)
    dirs = re.split(r'[/\\]', splitext[0])
    newdir = []
    for level in dirs:
        newdir.append(legalize(level))
    newdir[-1] = newdir[-1] + splitext[1]
    return os.path.join(*newdir)


def legalize(string: str):
    return re.sub(r'[^A-Za-z0-9_]', '_', string)


def dir_with_time(prefix: str):
    return prefix + '_' + str(time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()))
