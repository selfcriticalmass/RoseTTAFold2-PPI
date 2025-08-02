import os, sys

inlist = sys.argv[1] # a file containing protein pairs to be analyzed, e.g., protein_pairs
indir = sys.argv[2] # input directory containing single protein MSAs, e.g., protein_single_msas/
outdir = sys.argv[3] # output directory to store paired MSAs, e.g., protein_paired_msas/

fp = open(inlist,'r') 
pairs = []
for line in fp:
    words = line.split()
    pairs.append(words[0])
fp.close()

logp = open(inlist + '_input','w')
for pair in pairs:
    prot1 = pair.split('_')[0]
    prot2 = pair.split('_')[1]

    fp = open(indir + '/' + prot1 + '.a3m','r')
    acc2seqA = {}
    accsA = set([])
    for line in fp:
        if line[0] == '>':
            name = line[1:].split()[0]
        elif name == 'mask':
            maskA = line[:-1]
        elif name == 'query':
            qseqA = line[:-1]
        else:
            acc2seqA[name] = line[:-1]
            accsA.add(name)
    fp.close()

    good_indsA = set([])
    for i, char in enumerate(maskA):
        if char == '*':
            good_indsA.add(i)
    filt_qseqA = ''
    for i, char in enumerate(qseqA):
        if i in good_indsA:
            filt_qseqA += char

    filt_acc2seqA = {}
    for acc in accsA:
        seqA = acc2seqA[acc]
        ind = -1
        filt_seqA = ''
        for char in seqA:
            if char.islower():
                if ind in good_indsA and ind + 1 in good_indsA:
                    filt_seqA += char
            else:
                ind += 1
                if ind in good_indsA:
                    filt_seqA += char
        filt_acc2seqA[acc] = filt_seqA

    fp = open(indir + '/' + prot2 + '.a3m','r')
    acc2seqB = {}
    accsB = set([])
    for line in fp:
        if line[0] == '>':
            name = line[1:].split()[0]
        elif name == 'mask':
            maskB = line[:-1]
        elif name == 'query':
            qseqB = line[:-1]
        else:
            acc2seqB[name] = line[:-1]
            accsB.add(name)
    fp.close()

    good_indsB = set([])
    for i, char in enumerate(maskB):
        if char == '*':
            good_indsB.add(i)
    filt_qseqB = ''
    for i, char in enumerate(qseqB):
        if i in good_indsB:
            filt_qseqB += char

    filt_acc2seqB = {}
    for acc in accsB:
        seqB = acc2seqB[acc]
        ind = -1
        filt_seqB = ''
        for char in seqB:
            if char.islower():
                if ind in good_indsB and ind + 1 in good_indsB:
                    filt_seqB += char
            else:
                ind += 1
                if ind in good_indsB:
                    filt_seqB += char
        filt_acc2seqB[acc] = filt_seqB

    rp = open(outdir + '/' + pair + '.a3m','w')
    rp.write('>query\n')
    rp.write(filt_qseqA + filt_qseqB + '\n')
    for acc in accsA.intersection(accsB):
        seqA = filt_acc2seqA[acc]
        seqB = filt_acc2seqB[acc]
        rp.write('>' + acc + '\n')
        rp.write(seqA + seqB + '\n')
    rp.close()
    os.system('hhfilter -i ' + outdir + '/' + pair + '.a3m -o ' + outdir + '/' + pair + '.i90.a3m -id 90 >& tmplog')
    os.system('rm ' + outdir + '/' + pair + '.a3m')
    logp.write(outdir + '/' + pair + '.i90.a3m\t' + str(len(filt_qseqA)) + '\n')
logp.close()
