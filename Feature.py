from Bio.GenBank import Record


class Feature(Record.Feature):
    """
    child class of Genbank.Record.Feature

        Attributes:
         - feature - instance of parent feature class

        """

    def __init__(self, feature: Record.Feature):
        super().__init__()
        self.__dict__ = feature.__dict__

        # BioPython定义了一个qualifier类，里面含有两个字符串变量：key和value，然后把一堆qualifier类储存在一个list里面。
        # 而且，key和value里面还包括了原始文件的分隔符，例如，原始文件的原文为：/gene="trnM-CAU"
        # 储存为：key='/gene=', value='"trnM-CAU"'(斜杠、等于号、双引号都是字符串内的一部分
        # 所以我这里将其储存在一个qualifier_dict的字典里面，并且去掉斜杠、等于号、双引号等累赘的分隔符
        qualifier_dict = {}
        for qualifier in self.qualifiers:
            qualifier_dict[qualifier.key[1:-1]] = qualifier.value[1:-1]
        self.qualifier_dict = qualifier_dict

        # 对intervals或者叫location的储存是直接将原始文件的原文储存为一个字符串。
        # 例如：join(complement(70173..70286),140516..140746,complement(98372..98398))
        # 还有complement在外，join在内的：complement(join(103716..103750,104294..104525,76550..76663))
        # 此处调用后面构造的Interval类，将每一段作为一个interval，存储在一个list里面.
        # Interval类详见后文
        interval_strs = []
        self.intervals = []

        # 去掉外层的join和order，虽然我没搞明白这俩的区别是什么，但是……就这样吧。
        fully_complement = False
        if self.location.startswith("complement"):
            self.location = self.location[11:-1]
            fully_complement = True

        if self.location.startswith("join"):
            interval_strs.extend(self.location[5:-1].split(','))
        elif self.location.startswith("order"):
            interval_strs.extend(self.location[6:-1].split(','))
        else:
            interval_strs.append(self.location)

        if fully_complement:
            interval_strs.reverse()

        for interval_str in interval_strs:
            interval_str = interval_str.replace('<', "").replace(">", "")
            complement: bool
            if interval_str.startswith("complement"):
                # 如果以complement开头，去掉外层的complement，并将该段interval标记为反向互补序列
                interval_str = interval_str[11:-1]
                complement = not fully_complement  # 正常情况下为True，fully_complement情况下为False
            else:
                complement = fully_complement  # 正常情况下为False，fully_complement情况下为True
            interval_list = interval_str.split('..')
            interval = Interval(interval_list[0], interval_list[1], complement)
            self.intervals.append(interval)


class Interval(object):
    """
    详见上文，此处只考虑位置前后，位置在前者为start，位置在后者为end，不会因反向而颠倒。
    start和end是从1开始数的，所以对sequence取子串时需将start-1，
    而此处的end是指子串的最后一个字符，本来也需要-1，但是取子串的时候的end是子串结束后的下一个字符，需要+1，刚好抵消了，所以无需额外操作。
    例如：record.sequence[interval.start - 1:interval.end]

        Attributes:
         - start - 序列起始位置 start position
         - end - 序列终止位置 end position
         - complement - 是否反向互补序列 revert complement or not

        """

    def __init__(self, start: int, end: int, complement: bool):
        self.start = int(start)
        self.end = int(end)
        self.complement = complement

    def __str__(self):
        if self.complement:
            return '(' + str(self.end) + '..' + str(self.start) + ')'
        else:
            return '(' + str(self.start) + '..' + str(self.end) + ')'
