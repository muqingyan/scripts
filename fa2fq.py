import argparse

def run(args):
    qual = args.quality_score
    f = open(args.input)
    fout = open(args.output,'w')
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq != '':
                fout.write(seq+'\n')
                fout.write('+\n')
                fout.write(qual*len(seq)+'\n')
            fout.write('@'+line[1:])
            seq = ''
        else:
            seq += line.strip()
    f.close()
    fout.write(seq+'\n')
    fout.write('+\n')
    fout.write(qual*len(seq)+'\n')

def main():
    parser = argparse.ArgumentParser(description="Convert a fasta file to a fastq File")
    parser.add_argument('-in',help="fasta input file",dest="input",type =str,required=True)
    parser.add_argument('-out',help="fastq output file",dest="output",type =str,required=True)
    parser.add_argument('-qual',help="Quality score to fill in",dest="quality_score",type =str,default='?')
#    parser.set_defaults(func=run)
    args = parser.parse_args()
    run(args)
#    args.func(args)

if __name__ == "__main__":
    main()
