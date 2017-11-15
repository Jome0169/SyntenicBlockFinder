import getopt
import sys
import os
import re
import copy
from ParseGff import *
from VisualizeSyntenyBlocks import * 
from operator import itemgetter
from datetime import datetime



def SaveColinFile(arg1):
    """ Reads in collineartiy file and saves output into dictionary where key
    is the alignment, and value is a nested list of all genes hit.

    :arg1: Colliniearity file 
    :returns: Nested Dictionary

    '## Alignment 744: score=264.0 e_value=2.1e-08 N=6 cg05&pi2 minus': [['744-
    0:', 'ClCG05G017130.1', 'evm.model.Chr2.2290', '  4e-84'], 
    ['744-  1:', 'ClCG05G017180.1', 'evm.model.Chr2.2288', ' 4e-129'], 
    ['744-2:', 'ClCG05G017290.1', 'evm.model.Chr2.2282', '  4e-83'],
    ['744-  3:', 'ClCG05G017470.1', 'evm.model.Chr2.2271', '  1e-50']
    """
    DictOfAlgns = {}
    with open(arg1, 'r') as f:
        for _ in range(11):
            next(f)
        for line in f:
            if line.startswith("## Alignment"):
                AlgnName = line.strip()
                DictOfAlgns[AlgnName] = []
            else:
                GeneInfo = line.strip().split('\t')
                DictOfAlgns[AlgnName].append(GeneInfo)
    return DictOfAlgns

def ReadBaseName(arg1):
    """Reads in the base name of the proteins of the query proteins. This is
    used for parsing.

    :arg1: File containing base names that look like :
    => ../../00.src/CM3.5.1.protein.fa <==
    MEL
    => ../../00.src/Cmaxima_v1.1.protein.fa <==

    :returns: List of Base Names 
    ['MEL', 'Cma', 'Cmo', 'Cuc', 'Lsi', 'ClC', 'Csa', 'Cla']

    """
    FileBaseNames = []
    with open(arg1,'r') as f:
        for line in f:
            if line.startswith('='):
                pass
            else:
                name = line.strip()
                FileBaseNames.append(name)
    return FileBaseNames


def CreateLogicalAlignment(RawColinDictionary, PreFlag):
    """TODO: Docstring for CreateLigcalAlignmen.

    :t(RawColinDictionary): Raw Collinear file dictioanry from Save colin file
    :returns: TODO

    """
    QueryAlignments = {}
    counter = 0
    
    setlist = [] 
     
    for alignmetn, genegroup in RawColinDictionary.items():
        aligkey = "alg"
        geneorder = set()
        for list1 in genegroup:
            if PreFlag in list1[1]:
                geneorder.add(list1[1])
            elif PreFlag in list1[2]:
                geneorder.add(list1[2])
        algnname = aligkey + str(counter)
        QueryAlignments[algnname] = geneorder
        counter += 1
        setlist.append(geneorder)

    return QueryAlignments



def CreateFastAlignment(RawColinDictionary, PreFlag):
    """TODO: Docstring for CreateLigcalAlignmen.

    :t(RawColinDictionary): Raw Collinear file dictioanry from Save colin file
    :returns: TODO

    """

    QueryGenehit = {}

    for alignmetn, genegroup in RawColinDictionary.items():
        geneorder = []
        for list1 in genegroup:
            if PreFlag in list1[1]:
                MicroSet = (list1[1], list1[2])
                geneorder.append(MicroSet)
            elif PreFlag in list1[2]:
                MicroSet = (list1[2], list1[1])
                geneorder.append(MicroSet)
        QueryGenehit[alignmetn] = geneorder
   
    return QueryGenehit

def CreateOutputDirectory():
    """TODO: Docstring for CreateOutputDirectory.
    :returns: TODO

    """
    if not os.path.exists("Output"):
        os.makedirs("Output")
    else:
        pass



def ReadInBlastRepeatFile(RepeatFileName):
    """TODO: Docstring for ReadInBlastRepeatFile.

    :RepeatFileName: TODO
    :returns: TODO

    """
    ReadRepeatFile = []
    
    DictLocations = {}

    with open(RepeatFileName, 'r') as f:
        for line in f:
            Clean = line.strip().split()
            GenomeLocation = Clean[1]
            if GenomeLocation not in DictLocations:
                DictLocations[GenomeLocation] = [Clean]
            elif GenomeLocation in DictLocations:
                DictLocations[GenomeLocation].append(Clean)
    return DictLocations


class MCScanAlign(object):

    """Docstring for MCScanAlign. """

    def __init__(self, alignmentNumber, genesetlist):
        """TODO: to be defined1.

        :alignmentNumber: TODO
        :genesetlist: TODO

        """
        self._alignmentNumber = alignmentNumber
        self._genesetlist = genesetlist
        self._AlignmentIntersections = None
        self._BaseCounterAlign = None
        self._AllSubgenomes = True
        self._SyntBlocks = None
        self._FewMissingGenes = None
        self._GenomicLocations = None
        self._GenomicSequences = None
        
    

    def SearchSpeedDictionary(self, SpeedDictAlign):
        """TODO: Docstring for SearchSpeedDictionary.
        :returns: TODO

        """
        
        CollectAlgnHits = [] 
        for algn, set1 in  SpeedDictAlign.items():
            AlgnSpecList = [[algn]]
            for item in set1:
                if item[0] in self._genesetlist:
                    AlgnSpecList.append(item)
                else:
                    pass
            if len(AlgnSpecList) > 1:
                CollectAlgnHits.append(AlgnSpecList)
            else:
                pass
        self._AlignmentIntersections = CollectAlgnHits
        return self._AlignmentIntersections
    
    def ParseAlignmensByType(self, BasicNames, BasicRefName):
        """lets see if we have intersections in all 
        :returns: TODO

        """
        CounterDict4BaseNames = {}
        for item in BasicNames:
            CounterDict4BaseNames[item] = 0

        for item in CounterDict4BaseNames:
            if item in BasicRefName:
                CounterDict4BaseNames[item] += 1

        #'## Alignment 19: score=689.0 e_value=2.1e-33 N=14 0&wm2 minus',
        #('MELO3C026826T1', 'Cla013369'),
        for genehit in self._AlignmentIntersections:
            TakeSecondName = genehit[1][1]
            BasName = TakeSecondName[0:3]
            CounterDict4BaseNames[BasName] += 1 
        
        self._BaseCounterAlign = CounterDict4BaseNames
        return self._BaseCounterAlign


    def CheckForUnivesalAligbment(self):
        """TODO: Docstring for CheckForUnivesalAligbment.

        :f: TODO
        :returns: TODO

        """
        
        AllSubGenomesPresent = True
        for key, value in self._BaseCounterAlign.items():
            if value == 0:
                AllSubGenomesPresent = False    
            else:
                pass

        self._AllSubgenomes = AllSubGenomesPresent
        return self._AllSubgenomes

    def FormatAlignmentAndMatches(self):
        """TODO: Docstring for FormatAlignmentAndMatches.
        :returns: TODO

        """
        PackageAllAlignments = []
        FinalSyntblock = []
        Headrs = [self._alignmentNumber]
        HeadernameStorage = []

        for gene in self._genesetlist:
            Matches = [gene,]
            Headrs = [self._alignmentNumber]
            Counter = 1
            for valy in self._AlignmentIntersections:

                VsAlgn = "AlgnVs" + str(Counter)
                #Idk if i will ever need it, but i liek have it
                HeadernameStorage.append([valy[0][0],VsAlgn])
                Headrs.append(VsAlgn)
                Counter += 1
                Matchedgene = None
                for pair in valy[1:]:
                    if pair[0] == gene:
                        Matchedgene = pair[1]
                    else:
                        pass
                
                
                if Matchedgene != None:
                    Matches.append(Matchedgene)
                elif Matchedgene == None:
                    Matches.append("||")

            FinalSyntblock.append(Matches)
         

        FinalSyntblock.insert(0, Headrs)
        self._SyntBlocks = FinalSyntblock
        return self._SyntBlocks 



    def FilterForBestAlgn(self):
        """filters out bad alignments from the isolated syntentic blocks. We
        are going to have many many alignments, so not all deserve to pass
        certainly.
        :returns: TODO

        """
        
        NumberOfGenesInAlgn = len(self._genesetlist)
        NumberOfAlignments = len(self._SyntBlocks[0]) - 1

        Counter = 0
        NumberOfLines = 0
        
        for GeneHits in self._SyntBlocks[1:]:
            NumberOfLines += 1 
            for hit in GeneHits:
                if hit == "||":
                    Counter += 1

                elif hit != '||':
                    pass
        CalculateCellNumber = NumberOfAlignments  * NumberOfGenesInAlgn
        CalculatePercentageMissing = (float(Counter) / float(CalculateCellNumber))

        if CalculatePercentageMissing < .25:
            self._FewMissingGenes = True
            return self._FewMissingGenes 
        
        else:
            self._FewMissingGenes = False
            return self._FewMissingGenes 
        
    


    def WriteToOutPut(self, AlignmentFile):
        """TODO: Docstring for WriteToOutPut.
        :returns: TODO

        """
        with open(AlignmentFile, 'a') as f:

            for Header in self._SyntBlocks[0]:
                f.write('{:^20}'.format(Header))
            f.write('\n')

            for item in self._SyntBlocks[1:]:
                for thing in item:
                    f.write('{:^21}'.format(thing))
                f.write('\n')
            f.write('\n')
        #for name in Headrs:
        #    print('{:^18}'.format(name), end='')
        #print ('\n')
        #for item in FinalSyntblock:
        #    for thing in item:
        #        print('{:^18}'.format(thing), end='')
        #    print('\n')
        #print('\n')


    def ExtractGenomicRegions(self, InputGeneDictionary):
        """TODO: Docstring for ExtractGenomicRegions.

        :InputGeneDictionary: TODO
        :returns: TODO

        """

        CreateDiffList = []  

        FinalLocationset = []
        DictKeyset = InputGeneDictionary.keys()
        
        LenToIterateThrough = range(0,len(self._SyntBlocks[1]))
        for iternumber in LenToIterateThrough:
            MicroList = [item[iternumber]for item in self._SyntBlocks[1:]]
            CreateDiffList.append(MicroList)

        for genelist in CreateDiffList:
            BaseName = next(i for i in genelist if i != '||' )
            BaseStorage = BaseName[0:3]
            for Keyset in DictKeyset:
                if BaseStorage in Keyset:
                    Z = InputGeneDictionary[Keyset]
                    GeneRefAndHist = [Keyset]
                    for item in genelist:
                        for Genes in Z[:-1]:
                            for gene in Genes:
                                if item in gene:
                                    GeneRefAndHist.append(gene)
                    FinalLocationset.append(GeneRefAndHist)

        
        #CheckLensHere at some point
            



        self._GenomicLocations = FinalLocationset
        return self._GenomicLocations


    def ExtractNucloTides(self, LoadedFiles):
        """TODO:THIS CAPTURES REDUNDENT NTs AND MUST BE FIXED****
        :returns: TODO

        """
        SeqList = []
        IterThroughLocations = self._GenomicLocations


        for item in IterThroughLocations:
            
            ScaffPosition = item[1][0]
            KeySet = item[0]

            AllStarts = [int(val[2]) for val in item[1:]] 
            AllEnds = [int(val[3]) for val in item[1:]] 
            
            #for thing in item:
            #    print (thing)
            
            StartSeq = (int(min(AllStarts)))
            EndSeq = (int(max(AllEnds)))
            

            SequenceAssociated = LoadedFiles[KeySet][-1][ScaffPosition][StartSeq:EndSeq]
            MicroList = [item[0],ScaffPosition,str(StartSeq),str(EndSeq),SequenceAssociated]
            SeqList.append(MicroList)
        
        self._GenomicSequences = SeqList
        return self._GenomicSequences

        



    def PrintFinalOutput(self, Output1, Output2):
        """TODO: Docstring for PrintFinalOutput.

        :arg1: TODO
        :returns: TODO

        """
        Chunker = "###################################"
        
        with open(Output1, 'a') as f:
            ForMatAlgn =  Chunker +  self._alignmentNumber + Chunker
            f.write(ForMatAlgn)
            f.write('\n')

            for item in self._GenomicLocations:
                FileFoundIn = item[0][0]
                for thing in item[1:]:
                    StringConversion = [ str(x) for x in thing ]
                    FormatLine = FileFoundIn + '\t' + \
                    '\t'.join(StringConversion)
                    f.write(FormatLine)
                    f.write('\n')

        with open(Output2, 'a') as z:
            ForMatAlgn =  Chunker +  self._alignmentNumber + Chunker
            z.write(ForMatAlgn)
            z.write('\n')
            for SeqInfo in self._GenomicSequences:
                Header = ">" + SeqInfo[0][0] + '_' +'_'.join(SeqInfo[1:4])
                z.write(Header)
                z.write('\n')
                z.write(str(SeqInfo[-1].seq))
                z.write('\n')



def usage():
    print ("Usage: python FindSharedColinChunk.py ")
    print ("This program is to parse pairwise comparison outputs created by MCscanX.")
    print ("this is to make it easer to find sytenic blocks between multiple species you may")
    print ("be interested in")


def main():

    Iflag = None
    BFlag = None
    CFlag = None
    PreFlag = None
    RFlag = None
    Oflag = None
    

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:c:b:p:r:o:h", ["input",
            "basenames", "collineage", "prefix","reapeat", "output", "help"])

    except getopt.GetoptError:
        usage()

        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "-input"):
            Iflag = arg
        elif opt in ("-b", '-basenames'):
            BFlag = arg
        elif opt in ("-c", '-collinear'):
            CFlag = arg
        elif opt in ("-p", '-prefix'):
            PreFlag = arg
        elif opt in ("-r", '-repeat'):
            RFlag = arg
        elif opt in ("-o", "-output"):
            Oflag = arg
        elif opt in  ("-h", "-help"):
            usage()
            sys.exit(2)
        else:
           print ("Unhandeled options %s" % (otps))

    if Iflag == None:
        print ("Need Input list of files")
        usage()
        sys.exit(2)
    elif Oflag == None:
        print ("Need output file base name to write to ")
        usage()
        sys.exit(2)
    elif BFlag == None:
        print ("Need base file names for compared individusl(")
        sys.exit(2)
    elif CFlag == None:
        print ("need collinear file from MXCscan(")
        sys.exit(2)

    StartTime = datetime.now()
    #Remove OutFile if it already exists

    SequenceFile = str(Oflag) + ".algn.fasta"
    SequenceLocation = str(Oflag) + ".gene.loc.txt"
    CollinFile = str(Oflag) + ".combined.txt"
    
    try:
        os.remove(SequenceFile)
        os.remove(CollinFile)
        os.remove(SequenceLocation)
    except OSError:
        pass

    #Function Calls
    SavedColinFile = SaveColinFile(CFlag)
    QueryGenomePreFix = ReadBaseName(BFlag)
    QueryAlignments = CreateLogicalAlignment(SavedColinFile, PreFlag)
    FastAlignDict = CreateFastAlignment(SavedColinFile, PreFlag)
    
    print("Loading in RepeatLocaiton blast")
    RepeatLocations = ReadInBlastRepeatFile(RFlag)
    #LoadFilesIn
    print("Loadinging In Gff and Genomes")
    AllLoadedFiles = ParseFastaAndGffFile(BFlag) 


    #Create Class for Each Alignment, and filer
    print("Running Through Gene Class")
    for algn, geneSet in QueryAlignments.items():
        AlgnObj = MCScanAlign(algn, geneSet)
        AlgnObj.SearchSpeedDictionary(FastAlignDict)
        AlgnObj.ParseAlignmensByType(QueryGenomePreFix,PreFlag)
        Z = AlgnObj.CheckForUnivesalAligbment()
        #If one hit from all sub genomes in algnt, proceed, else fail
        if Z == True:
            AlgnObj.FormatAlignmentAndMatches()
            AlgnObj.FilterForBestAlgn()
            if AlgnObj._FewMissingGenes == True:
                AlgnObj.ExtractGenomicRegions(AllLoadedFiles)
                AlgnObj.ExtractNucloTides(AllLoadedFiles)
                #AlgnObj.WriteToOutPut(CollinFile)
                #AlgnObj.PrintFinalOutput(SequenceLocation,SequenceFile)
                print("Entering Synteny Drawer")
                SyntenyDrawer(AlgnObj._GenomicLocations, AlgnObj._SyntBlocks,\
                        AlgnObj._alignmentNumber, RepeatLocations, \
                        AlgnObj._GenomicSequences)
            elif AlgnObj._FewMissingGenes == False:
                pass

        elif Z == False:
            pass

       
    #Speed Things
    EndTime = datetime.now()
    FinalTime = EndTime - StartTime

    print ("Total Time %s" % (FinalTime))

if __name__=="__main__":
    main()
