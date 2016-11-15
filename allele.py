class Allele(object):
    """Abstraction of an IGHV allele
        This class will accumulate occurrences, D and J assignments 

    Attributes:
        name : name of this allele
        vgene: IGHV gene assignment
        dass : IGHD gene assignments
        jass : IGHJ gene assignments
        cdr3s: set of CDR3s associated with this allele 
        nseqs: # of sequences assigned to this allele
        exact: # of exact matches
    """

    def __init__(self, name):
        self.name  = name
        self.vgene = name.split('*')[0]
        self.dass  = set()
        self.jass  = set()
        self.cdr3s = set()
        self.nseqs = 0
        self.exact = 0

    def addRow(self, count, d, j, cdr3aa, exact_bool):
        """ Adds a row of IgBLAST hits into this allele """

        self.dass.add(d)
        self.jass.add(j)
        self.cdr3s.add(cdr3aa)

        self.nseqs += int(count)
        if exact_bool:
            self.exact += int(count)

    def shortname(self):
        return self.name.split('*')[1]

    def snr(self):
        return shortname().split('_')[1]

    def __eq__(self,other):
        return self.name == other.name

    def __str__(self):
        return "{0} allele of {1} gene".format(self.shortname(), self.vgene)
    
    def __repr__(self):
        return "{0} [# of seq = {1}, Exact = {2}, Ds = {3}, Js = {4}, CDR3s = {5}]" \
            .format(self.name, self.nseqs, self.exact, len(self.dass), len(self.jass), len(self.cdr3s))
       
    def __hash__(self):
        return hash(self.name)
