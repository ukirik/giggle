import csv

class STranslator(object):
    """ A utility class to translate the S-codes to IMGT names"""

    def __init__(self, reffile):
        self.d = {}
        with open(reffile) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.d[row["IgDiscover_S-kod"]] = row["Germline"].strip()

    def get(self, scode):
        if scode not in self.d:
            print("Warning, code {0} not found!".format(scode))
            return scode
        return self.d[scode]

