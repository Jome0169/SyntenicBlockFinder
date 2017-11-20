from sys import argv
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from operator import itemgetter
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
import random



def SyntenyDrawer(LocationData, SytenyBlockInfo,alignmentname, \
        RepeatLocations, GenomeSeqLocations):
    """TODO: Docstring for SyntenyDrawer.
    :returns: TODO

    """
        
    def SyntBlockColr(SyntInfo):
        """Takes in syntenic block in formation and returns the colors to be
        associated with each block, as well as dictionary associated with the
        cross link between genes. Generates based off number of input genes and
        randomly generates the colors.

        :SyntInfo: TODO
        :returns: TODO

        """
        GenomeColors = {}
        ColorTown  = []

        ColorCount = 100
        while ColorCount != 0:
            X = "#%06x" % random.randint(0, 0xFFFFFF)
            ColorTown.append(colors.HexColor(X))
            ColorCount -= 1 

        LenToIterateThrough = range(0,len(SyntInfo[0]))
        Counter = 0
        for item in SyntInfo[1:]:
            GenomeColors[ColorTown[Counter]] = item
            Counter += 1 
        return GenomeColors 
    
    
    def QuickParseRepeats(genomiclocationstring, TEDict, NucStartLoc):
        TEPairs = []
        

        Header = genomiclocationstring[0][0] + '_' \
        +'_'.join(genomiclocationstring[1:4])

        z = list(TEDict.keys())
        if Header in TEDict:
            for item in TEDict[Header]:
                if int(item[8]) < int(item[9]):
                    TakePosition = [int(item[8]), int(item[9])]
                    TEPairs.append(TakePosition)
                elif int(item[8]) > int(item[9]):
                    TakePosition = [int(item[9]), int(item[8])]
                    TEPairs.append(TakePosition)
        return TEPairs


            #for SeqInfo in self._GenomicSequences:
            #    Header = ">" + SeqInfo[0][0] + '_' +'_'.join(SeqInfo[1:4])


    X = SyntBlockColr(SytenyBlockInfo)
    Largest = 0 

    gdd = GenomeDiagram.Diagram('Diagram Name')
    #Run through at the same time to ADD TEs based off Header Location info
    for item, seqinfo in zip(LocationData, GenomeSeqLocations) :
        
        #Find All Startings and endings
        AllStarts = [int(val[2]) for val in item[1:]] 
        AllEnds = [int(val[3]) for val in item[1:]] 
        StartSeq = (int(min(AllStarts)))
        EndSeq = (int(max(AllEnds)))
    

        #FindTEs
        print("FInding TE Pairs")
        FoundTEPairs = QuickParseRepeats(seqinfo,RepeatLocations, StartSeq)
        

        print("Drawing in Progress")
        #Initialize features 
        gdt_features = gdd.new_track(1, greytrack=True, name=[item[0][0]],
                start=0, end=EndSeq)
        gds_features = gdt_features.new_set()
      
        CopyList = item[1:]
        for list1 in CopyList:
            ReformatStart = int(list1[2]) - StartSeq
            ReformatEnd = int(list1[3]) - StartSeq 
            list1[2] = ReformatStart
            list1[3] = ReformatEnd

            if ReformatEnd > Largest:
                Largest = ReformatEnd 
            else:
                pass

        sorted(CopyList, key=itemgetter(2))
        for NewFeat in CopyList:
            for key,val in X.items():
                if NewFeat[1] in val:
                    Q = key

            feature = SeqFeature(FeatureLocation(NewFeat[2], \
            NewFeat[3]), strand=None)
            gds_features.add_feature(feature, color=Q, name=str(NewFeat[1]), \
                    label=True,label_position="middle", label_size = 6,\
                    label_angle=-90 )

        if len(FoundTEPairs) != 0:
            print(item[0], alignmentname, len(FoundTEPairs))
            for TE in FoundTEPairs:
                TEStarts = TE[0] 
                TEEnd = TE[1] 
                print(TE, TEStarts, TEEnd)
                print(CopyList)
                
                feature = SeqFeature(FeatureLocation(TEStarts, TEEnd), strand=None)
                gds_features.add_feature(feature,
                        color="purple", name=str("MITE"), label=True, \
                        label_size = 6, label_angle=0)
        else:
            pass


    NameDraw = str(alignmentname) + "_genomediagramtest.pdf"
    gdd.draw(format='linear', pagesize=(30*cm,30*cm), fragments=1,
                         start=0, end=Largest)
    gdd.write(NameDraw, "pdf")
