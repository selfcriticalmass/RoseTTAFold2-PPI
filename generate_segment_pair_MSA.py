import os, sys

inlist = sys.argv[1] # a file containing segment pairs to be analyzed, e.g., segment_pairs
indir = sys.argv[2] # input directory containing single segment MSAs, e.g., segment_single_msas/
outdir = sys.argv[3] # output directory to store paired MSAs, e.g., segment_paired_msas/

fp = open(inlist,'r')
pairs = []
for line in fp:
    words = line.split()
    pairs.append(words[0])
fp.close()

logp = open(inlist + '_input','w')
for pair in pairs:
    seg1 = pair.split('__')[0]
    seg2 = pair.split('__')[1]

    fp = open(indir + '/' + seg1 + '.a3m','r')
    acc2seqA = {}
    accsA = set([])
    for line in fp:
        if line[0] == '>':
            name = line[1:].split()[0]
        elif name == 'query':
            qseqA = line[:-1]
        else:
            acc2seqA[name] = line[:-1]
            accsA.add(name)
    fp.close()

    fp = open(indir + '/' + seg2 + '.a3m','r')
    acc2seqB = {}
    accsB = set([])
    for line in fp:
        if line[0] == '>':
            name = line[1:].split()[0]
        elif name == 'query':
            qseqB = line[:-1]
        else:
            acc2seqB[name] = line[:-1]
            accsB.add(name)
    fp.close()

    rp = open(outdir + '/' + pair + '.a3m','w')
    rp.write('>query\n')
    rp.write(qseqA + qseqB + '\n')
    for acc in accsA.intersection(accsB):
        seqA = acc2seqA[acc]
        seqB = acc2seqB[acc]
        rp.write('>' + acc + '\n')
        rp.write(seqA + seqB + '\n')
    rp.close()
    os.system('hhfilter -i ' + outdir + '/' + pair + '.a3m -o ' + outdir + '/' + pair + '.i90.a3m -id 90 >& tmplog')
    os.system('rm ' + outdir +'/' + pair + '.a3m')
    logp.write(outdir + '/' + pair + '.i90.a3m\t' + str(len(qseqA)) + '\n')
logp.close()
