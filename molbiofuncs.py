''' this is a module of utility functions for routine computational molecular biology manipulations '''

geneticcode = {'ggg':'G', 'gga':'G', 'ggc':'G', 'ggt':'G','gag':'E', 'gaa':'E', 'gac':'D', 'gat':'D','gcg':'A', 'gca':'A', 'gcc':'A', 'gct':'A',
               'gtg':'V', 'gta':'V', 'gtc':'V', 'gtt':'V','agg':'R', 'aga':'R', 'agc':'S', 'agt':'S','aag':'K', 'aaa':'K', 'aac':'N', 'aat':'N',
               'acg':'T', 'aca':'T', 'acc':'T', 'act':'T','atg':'M', 'ata':'I', 'atc':'I', 'att':'I','cgg':'R', 'cga':'R', 'cgc':'R', 'cgt':'R',
               'cag':'Q', 'caa':'Q', 'cac':'H', 'cat':'H','ccg':'P', 'cca':'P', 'ccc':'P', 'cct':'P','ctg':'L', 'cta':'L', 'ctc':'L', 'ctt':'L',
               'tgg':'W', 'tga':'*', 'tgc':'C', 'tgt':'C','tag':'*', 'taa':'*', 'tac':'Y', 'tat':'Y','tcg':'S', 'tca':'S', 'tcc':'S', 'tct':'S','ttg':'L', 'tta':'L', 'ttc':'F', 'ttt':'F'}

def translate(seq):
    '''converts DNA to Amino Acid single letter code\n
    takes an in-frame DNA sequence as the argument'''
    i = 0
    aa = ""
    while i <= len(seq)-3:
        try:
            aa += geneticcode[seq[i:i+3].lower()]
        except KeyError:
            aa += "X"
        i +=3
    return aa

def find_stop_codons(sequence):
    '''locates stop codons and returns them and thier location (1 based)\n
       returns a dictionary of codon:position pairs\n
       takes an in-frame DNA sequence as an argument'''
    i = 0
    stops = []
    locations = []
    while i <= len(sequence)-3:
        try:
            if "*" == geneticcode[sequence[i:i+3].lower()]:
                stops += [sequence[i:i+3]]
                locations.append(i+1)
        except KeyError:
            pass
        i +=3
    return dict(zip(stops, locations))

def rev_comp(seq):
    '''reverse complements a sequence'''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def rev(sequence):
    '''reverses a sequence'''
    return sequence[::-1]

def comp(seq):
    '''complements a sequence'''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
    return "".join(complement.get(base, base) for base in seq)

def getKmers(sequence, size=6, step=1):
    '''returns all possible k-mers for the given sequence\n
    takes a sequence and a k-mer size as arguments\n
    the default size = 6\n
    converts to all lower case\n
    specify the stride length (stagger between k-mers) with step argument'''
    return [sequence[x:x+size].lower() for x in range(0, len(sequence) - size + 1, step)]

#seq = 'ACGTTCAAGTACGATTAGGATCCACGAACCTGA'
#translate(seq)
#find_stop_codons(seq)
#rev_comp(seq)
#rev(seq)
#comp(seq)
#rev(comp(seq))
#getKmers(seq, size=6, step=2)

