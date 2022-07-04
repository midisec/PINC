from functools import reduce
import pandas as pd
import os
import sys, getopt

from autogluon.tabular import TabularPredictor
import pandas as pd
from sklearn.metrics import f1_score
from sklearn import metrics
from sklearn.metrics import roc_curve, auc


class Tools(object):

    def __init__(self, path1, feather, index_col, path2, result_path):
        self.path1 = path1
        self.feather = feather
        self.index_col = index_col
        self.path2 = path2
        self.mySeries = self.init1()
        self.result_path = result_path

    def init1(self):
        def concate_fullconnect(mer_set, bases):
            res = []
            for i in mer_set:
                for j in bases:
                    res.append(i + j)

            return res

        mylist = []
        mylist.append('seqLen')
        mylist.append('GC_C')
        mylist.append('score')
        mylist.append('cdsStarts')
        mylist.append('cdsStop')
        mylist.append('cdsSizes')
        mylist.append('cdsPercent')
        mylist.append('label')
        mylist3 = reduce(concate_fullconnect, [['A', 'T', 'C', 'G']] * 3)
        mylist2 = reduce(concate_fullconnect, [['A', 'T', 'C', 'G']] * 2)
        mylist1 = reduce(concate_fullconnect, [['A', 'T', 'C', 'G']] * 1)
        mylist = mylist + mylist1 + mylist2 + mylist3
        return pd.Series(mylist)

    def readFa(self, fa):
        '''
        @msg: Read a FASTA file
        @param fa {str}  fasta File path
        @return: {generator} Return a generator, which can iteratively get each sequence name and sequence of the FASTA file
        '''
        with open(fa, 'r') as FA:
            seqName, seq = '', ''
            while 1:
                line = FA.readline()
                line = line.strip('\n')
                if (line.startswith('>') or not line) and seqName:
                    yield ((seqName, seq))
                if line.startswith('>'):
                    seqName = line[1:].split()[0]
                    seq = ''
                else:
                    seq += line
                seq = seq.upper()
                if not line: break

    # Get sequence
    def getSeq(self, fa, querySeqName, start=1, end=0):
        '''
        @msg: Get a sequence of FASTA files
        @param fa {str}  FASTA file path
        @param querySeqName {str}  Sequence name
        @param start {int}  When intercepting this sequence, the starting position can be omitted, and the default is 1
        @param end {int}  When FASTA intercepts the sequence, the last position can be omitted, and the default is the full length of the sequence
        @return: {str} Return the sequence found (intercepted)
        '''
        if start < 0: start = start + 1
        for seqName, seq in readFa(fa):
            if querySeqName == seqName:
                if end != 0:
                    returnSeq = seq[start - 1:end];
                    print(start - 1)
                else:
                    returnSeq = seq[start - 1:]
                return returnSeq

    # Get reverse complementary sequence
    def getReverseComplement(self, sequence):
        '''
        @msg: Get reverse complementary sequence
        @param sequence {str}
        @return: {str} Return reverse complementary sequence
        '''
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()[::-1]

    # Get GC content
    def getGC(self, sequence):
        '''
        @msg: Get the GC content of a sequence
        @param sequence {str}
        @return: {float} Return GC content
        '''
        sequence = sequence.upper()
        return (sequence.count("G") + sequence.count("C")) / len(sequence)

    # Get window sequence
    def readSeqByWindow(self, sequence, winSize, stepSize):
        '''
        @msg: Sliding window reads a sequence
        @param sequence {str}
        @param winSize {int} Window size
        @param stepSize {int}  step
        @return: {generator}  A generator is returned, and each window sequence of the sequence can be iteratively obtained
        '''
        if stepSize <= 0: return False
        now = 0
        seqLen = len(sequence)
        while (now + winSize - stepSize < seqLen):
            yield sequence[now:now + winSize]
            now += stepSize

    # Get gap location
    def getGapPos(self, sequence):
        '''
        @msg: Get the position of gap in a sequence
        @param sequence {str}
        @return: {list}  Return a list. Each element in the list is the start and end position of each gap
        '''
        Ns = {'N', 'n'}
        result = []
        i = 0
        for base in sequence:
            i += 1
            if not base in Ns: continue
            if len(result) == 0:
                result.append([i, i])
            elif i - result[-1][1] == 1:
                result[-1][1] = i
            else:
                result.append([i, i])
        return result

    # Get k-mer
    def printSeq(self, sequence):
        # Extract your code into a function and print header for current kmer
        # print("%s" %name)
        dis = {}
        kmers = {}
        for k in range(1, 4):
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                if kmer in kmers:
                    kmers[kmer] += 1
                else:
                    kmers[kmer] = 1

        for kmer, count in kmers.items():
            # print (kmer + "\t" + str(count))
            dis[kmer] = str(count)
        return dis

    def task(self):
        # Get the text of CDs prediction tool (mRNA)
        rna_cds = pd.read_csv(self.path1,
                              names=self.feather,
                              sep='\t', index_col=self.index_col)

        seqName_list = []

        for seqName, seq in self.readFa(self.path2):
            seqLen = len(seq)
            seqName_list.append(seqName)

        Rna_f = pd.DataFrame(0, index=seqName_list, columns=self.mySeries)
        for seqName, seq in self.readFa(self.path2):
            seqLen = len(seq)
            seqName_list.append(seqName)
            GC = self.getGC(seq)
            dis = self.printSeq(seq)
            Rna_f.loc[seqName, 'GC_C'] = GC
            Rna_f.loc[seqName, 'seqLen'] = seqLen
            for k_mer in dis:
                Rna_f.loc[seqName, k_mer] = dis[k_mer]

        for index, row in rna_cds.iterrows():
            # test_1.loc[row["seqName"],'cdsStarts']=row["cdsStarts"]
            # print(index)
            # row.info
            tem_pd_score = pd.DataFrame(row)
            Rna_f.loc[index, 'cdsStarts'] = row["start"]
            Rna_f.loc[index, 'cdsStop'] = row["end"]
            Rna_f.loc[index, 'cdsSizes'] = row["end"] - row["start"]
            Rna_f.loc[index, 'score'] = tem_pd_score.iloc[4, 0]
            Rna_f.loc[index, 'cdsPercent'] = (row["start"] + row["end"]) / Rna_f.loc[index, 'seqLen']
            # test_1.loc[row["seqName"],"cdsStarts"]=row["cdsStarts"]

        a = [1, 2, 4, 5, 9, 10, 25, 26, 39, 45]

        #    todo file save location
        Rna_f.iloc[:, a].to_csv(self.result_path)

    def run(self):
        self.task()


class PreTools(object):
    def __init__(self, source_name):
        # self.source_path = "/app/pinc/data/sample_old/" + source_name + ".fasta"
        self.source_path = "./" + source_name + ".fasta"

        self.result_path = "/app/pinc/data/sample_pre/"
        self.result_name = source_name + ".cds"

    def check_path(self):
        pass

    def produce(self):
        cmdline = "/app/kent/kentUtils-master/bin/txCdsPredict {} {}".format(self.source_path,
                                                                             self.result_path + self.result_name)

        os.system(cmdline)

    def run(self):
        self.produce()


class Prediction(object):
    def __init__(self, result_path, name):
        self.path = "model"
        self.result_path = result_path
        self.name = name
        self.predictor = TabularPredictor.load(self.path)

        self.data = pd.read_csv(result_path)
        self.test = self.data
        # self.test = self.test.drop(columns=['label'])

    def run(self):
        proba = self.predictor.predict_proba(self.test)
        proba0 = proba[0]
        proba1 = proba[1]

        with open('/app/pinc/prediction_result/'+self.name+'.csv', 'w') as f:
            f.write("id,Label,Coding probability / 0,Non-coding probability / 1")
            f.write("\n")
            for i in range(len(self.test.iloc[:, 0])):
                # print("{},{},{},{}".format(self.test.iloc[:, 0][i], "Non-coding" if proba0[i] > proba1[i] else "coding", proba0[i], proba1[i]))
                f.write("{},{},{},{}".format(self.test.iloc[:, 0][i], "coding" if proba0[i] > proba1[i] else "Non-coding", proba0[i], proba1[i]))
                f.write("\n")

        count = 0
        for i in range(len(self.test.iloc[:, 0])):
            count = count + 1
            if count>50:
                print("[!] Only print the first 50 data items")
                break
            print("{},{},{},{}".format(self.test.iloc[:, 0][i], "Non-coding" if proba0[i] > proba1[i] else "coding", proba0[i], proba1[i]))


        print("[*] The result save path: {}".format('/app/pinc/prediction_result/'+self.name+'.csv'))


class OneForAll(object):
    def __init__(self, source_name):
        self.source_path = source_name

        self.names = ['seqname', 'start', 'end', 'source ', 'accession', 'score ', 'startComplete ',
                                 'endComplete ', 'cdsCount ', 'cdsStarts ', 'cdsSizes ']
        self.index_col = "seqname"

    def run(self):
        pre = PreTools(self.source_path)
        pre.run()

        t1 = Tools("/app/pinc/data/sample_pre/{}".format(self.source_path+".cds"), self.names, self.index_col, "./" + self.source_path + ".fasta", "/app/pinc/data/result/{}".format(self.source_path+".csv"))
        t1.run()

        prediction = Prediction("/app/pinc/data/result/{}".format(self.source_path+".csv"), self.source_path)
        prediction.run()


def print_logo():
    print("""
      ____ ___ _   _  ____ 
     |  _ \_ _| \ | |/ ___|
     | |_) | ||  \| | |    
     |  __/| || |\  | |___ 
     |_|  |___|_| \_|\____|                   
    """)
    print("usage: \n"
          "     pinc.py -h                  # View Help\n"
          "     pinc.py -f <data.fasta>     # Data Prediction\n")


if __name__ == '__main__':
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hf:c:", ["ifile=", "ofile="])
        if len(argv) == 0:
            print_logo()
            sys.exit(2)
    except getopt.GetoptError:
        print_logo()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_logo()
            sys.exit()
        elif opt in ("-f", "--ifile"):
            if "/" in arg:
                arg = arg.split("/")[-1]

            if ".fasta" != arg[-6:]:
                print_logo()
                print("ERROR: Incorrect file extensions")
                sys.exit(2)
            source_path = arg.split('.')[0]

    task = OneForAll(source_path)
    task.run()
