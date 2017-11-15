import sys
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os
from sys import argv






def ParseFastaAndGffFile(BaseFileNameSystem):
    """TODO: Docstring for ParseFastaAndGffFile.

    :BaseFileNameSystem: TODO
    :returns: TODO

    """

    def ReadGffFiles(FileName, DictionaryOfGffs):
        """TODO: Docstring for ReadGffFiles.

        :arg1: TODO
        :returns: TODO

        """
        GeneList = []
        FilteredGeneList = []
        
        GffFile = FileName[0] + ".gff3"
        with open(GffFile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    pass
                elif len(line) == 0:
                   pass 
                else:
                    Cleanedline = line.strip().split('\t')
                    if Cleanedline[0] == '':
                        pass
                    elif Cleanedline[2] == "mRNA":
                        GeneList.append(Cleanedline)
                    else:
                        pass

        for item in GeneList:
            ChrName = item[0]
            GeneName = item[8]
            if ':' in GeneName:
                GeneName = item[8].split(":")
            elif ';' in GeneName:
                GeneName = item[8].split(";")
            
            GeneNameIso = GeneName[0].replace('ID=','')
            GeneStart = str(int(item[3]) -1)
            GeneEnd = str(int(item[4]))
            Combined = [ChrName,GeneNameIso,GeneStart,GeneEnd]
            FilteredGeneList.append(Combined)

        DictionaryOfGffs[FileName] = [FilteredGeneList]
        

                 

    def ReadFastFile(FileName, DictionaryOfFiles):
        """TODO: Docstring for ReadFastFile.

        :arg1: TODO
        :returns: TODO

        """
        FastaFile = FileName[0] + '.genome.fa' 
        ChromScafParse = SeqIO.to_dict(SeqIO.parse(FastaFile, 'fasta'))
        DictionaryOfFiles[FileName].append(ChromScafParse)

    FileNameList = []
    FinalFileDict = {}
    BasPair = []

    with open(BaseFileNameSystem, 'r')  as f:
        BasPair = []
        for line in f:
            if line.startswith('=>'):
                CleanedLine = line.split(' ')
                X = [CleanedLine[1]]
                BaseNameOfFile = CleanedLine[1]
                FileNameList.append(BaseNameOfFile)
            else:
                CleanedLine = line.strip()
                X.append(CleanedLine)
                BasPair.append(X)


                

    for BaseName in BasPair:
        ReadGffFiles(tuple(BaseName), FinalFileDict)
        ReadFastFile(tuple(BaseName), FinalFileDict)

    return FinalFileDict









