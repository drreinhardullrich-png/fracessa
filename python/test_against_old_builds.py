import os
import multiprocessing
from pprint import pprint
# import functools
import time
# import re
import sys


# format: dimension, number of ess, is circular symmetric, matrix input
inputmatrices = [
    [2, 1, False, "0,1,1,0"],
    [2, 2, False, "3,1,2,4"],
    [3, 2, False, "4,6,1,7,5,4,0,7,3"],
    [3, 1,  False, "1,2,2,2,1,2,0,0,1"],
    [3, 1, False, "0,0,0,0,0,0,1,1,0"],
    [3, 1, False, "0,2,-1,2,0,3,-1,3,1/2"],
    [4, 0, False, "-1,1,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0"],
    [5, 5, True, "1,3"],
    [5, 6, False, "1,0,2,2,2,0,1,2,2,2,2,2,1,0,0,2,2,0,1,0,2,2,0,0,1"],
    [6, 9, True, "1,-1,1"],
    [7, 14, True, "1,5,8"],
    [7, 14, True, "8,13,2"],
    [8, 20, True, "5,5,13,11"],
    [9, 30, True, "1,5,10,13"],
    [10, 50, True, "13,19,25,7,5"],
    [11, 66, True, "1,4,6,6,5"],
    [11, 55, True, "15,1/100,14,11/10,76/10"],
    [12, 105, True, "533,118,326,357,118,477"],
    [12, 12, True, "287,50,200,150,113,200"],
    [12, 24, True, "292,-1,293,-1,293,-1"],
    [12, 105, True, "16029525,4007380,9088537,12022143,2147548,16029525"],
    [13, 143, True, "7,18,18,10,7,10"],
    [14, 224, True, "4,3,2,2,3,4,-1"],
    [15, 360, True, "9,9,5,9,3,5,9"],
    [16, 512, True, "7,5,4,4,4,5,7,-1"],
    [17, 493, True, "1627,1039,1527,1075,1527,1039,1627,648"],
    [18, 1152, True, "8,7,5,5,5,5,7,8,-1"],
    [19, 1444, True, "13,31,31,27,31,27,13,13,27"],
    [19, 19, True, "1,2,2,2,2,2,1,1,2"],
    [20, 2560, True, "13,13,13,8,8,8,13,13,13,-1"],
    [21, 4410, True, "15,15,7,15,15,7,7,15,7,15"],
    [22, 5632, True, "2,2,4,4,5,5,4,4,2,2,-1"],
    [23, 2507, True, "27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939"],
    [24, 15120, True, "15,15,7,15,15,7,15,7,7,15,15,7"]
]


class Tester:

    counter = 1

    def __init__(self, dimension, ess_expected, cs, matrix, typ, exe, parameters=''):
        self.id = Tester.counter
        Tester.counter += 1
        self.typ = typ
        self.exe = exe
        self.parameters = parameters
        self.dimension = dimension
        self.ess_expected = ess_expected
        self.cs = cs
        self.matrix = matrix
        self.name = '{}#{}'.format(self.dimension, self.matrix)
        self.syscall = '{} {}{} {}'.format(self.exe, '-' if self.parameters else '', self.parameters, self.name)
        self.result = 'xxx'
        self.header = ''
        self.vectors = list()

    def compute(self):
        # print(self.syscall)
        t0 = time.perf_counter()
        self.result = os.popen(self.syscall).read()
        t1 = time.perf_counter()
        self.timing = t1-t0
        # print(result)

        # data.result = re.sub(r';\d*\.*\d*\n', '\n', output)
        splitlines = self.result.splitlines()
        self.ess = int(splitlines[0])
        if len(splitlines) > 1:
            header = splitlines[1]
            vectors = splitlines[2:]
            header = header.split(';')
            # pprint(header)
            transformer = {
                'VectorID': 'candidate_id',
                'Vector': 'vector',
                'Support': 'support',
                'SupportSize': 'support_size',
                'ExtendedSupport': 'extended_support',
                'ExtendedSupportSize': 'extended_support_size',
                'ShiftReference': 'shift_reference',
                'IsEss': 'is_ess',
                'Reason': 'reason_ess',
                'Payoff': 'payoff',
                'PayoffDecimal': 'payoff_double'
            }

            for line in vectors:
                x = dict()
                for idx, item in enumerate(line.split(';')):
                    # if header[idx] not in ('PayoffDecimal', 'VectorID', 'ShiftReference'):
                    # x[header[idx]] = re.sub(r'/1$', '', item)
                    if self.typ == 'werner':
                        x[transformer[header[idx]]] = item
                    else:
                        x[header[idx]] = item
                self.vectors.append(x)

        print(self.syscall, self.timing)

    def check_against(self, other):

        if self.ess != other.ess:
            print(self.name, self.typ, other.typ, 'ess not equal: {} {}'.format(self.ess, other.ess))
            sys.exit()

        # pprint(self.vectors)
        if len(self.vectors) != len(other.vectors):
            print(self.name, self.typ, other.typ, 'len vectors not equal')
            pprint(self.__dict__)
            pprint(other.__dict__)

        used = set()
        for v in self.vectors:
            found = False
            for idx, ov in enumerate(other.vectors):
                if v == ov:
                    used.add(idx)
                    found = True
                    break
            if not found:
                print(self.name, self.typ, other.typ, 'vector not found: {}'.format(v))
                pprint(self.__dict__)
                pprint(other.__dict__)
                sys.exit()
        if len(used) != len(self.vectors):
            print('one vector was used double!')
            sys.exit


def input_testdata(iq, parameters_new, parameters_old, max_dimension):

    for line in inputmatrices:

        # x = Tester(line[0], line[1], line[2], line[3], 'mono /home/reinhard/REF/_EFR/EFR/EFR/bin/Release/EFR.exe', parameters)
        # x.typ = 'csharp'
        # testdata.append(x)

        x = Tester(line[0], line[1], line[2], line[3], 'werner', './REF_version_werner.exe', parameters_old)
        if x.dimension <= max_dimension:
            iq.put(x)

        x = Tester(line[0], line[1], line[2], line[3], 'github', '../build_linux64/REF_linux64.exe', parameters_new)
        if x.dimension <= max_dimension:
            iq.put(x)

    with open('testmatrizenfuerref_werner.txt') as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    for line in content:

        if not line:
            continue

        line = line.split('#')

        # x = Tester(int(line[0]), None, True, line[1], 'mono /home/reinhard/REF/_EFR/EFR/EFR/bin/Release/EFR.exe', parameters)
        # x.typ = 'csharp'
        # testdata.append(x)

        x = Tester(int(line[0]), None, True, line[1], 'werner', './REF_version_werner.exe', parameters_old)
        if x.dimension <= max_dimension:
            iq.put(x)

        x = Tester(int(line[0]), None, True, line[1], 'github', '../build_linux64/REF_linux64.exe', parameters_new)
        if x.dimension <= max_dimension:
            iq.put(x)


def worker(iq, oq):
    while True:
        x = iq.get()
        x.compute()
        oq.put(x)
        iq.task_done()


if __name__ == '__main__':

    parameters_new = 'c'
    parameters_old = 'v'
    max_dimension = 18

    iq = multiprocessing.JoinableQueue()
    oq = multiprocessing.Queue()

    input_testdata(iq, parameters_new, parameters_old, max_dimension)

    processes = []
    for i in range(multiprocessing.cpu_count()-2):
        worker_process = multiprocessing.Process(target=worker, args=(iq, oq), daemon=True, name='worker_process_{}'.format(i))
        worker_process.start()
        processes.append(worker_process)

    iq.join()

    testdata = list()
    while oq.qsize() != 0:
        testdata.append(oq.get())
    testdata = sorted(testdata, key=lambda k: k.id)

    print('############################# results ##################################')
    names = [x.name for x in testdata]
    names = list(set(names))
    for name in names:
        test = [x for x in testdata if x.name == name]
        first = test[0]
        for x in test[1:]:
            first.check_against(x)

    for x in testdata:
        print(x.typ, x.name, x.ess, x.timing)
        #     x.vectors = None
        #     pprint(x.__dict__)
        # print(x.name, len(x.vectors))
