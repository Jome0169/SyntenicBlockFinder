from sys import argv
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from operator import itemgetter
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
import random



def SyntenyDrawer(LocationData, SytenyBlockInfo,alignmentname):
    """TODO: Docstring for SyntenyDrawer.
    :returns: TODO

    """
        

    def SyntBlockColr(SyntInfo):
        """Takesi n syntenic block in formation and returns the colos to be
        associated with each block, as well as dictionary associated with the
        cross link between geenes

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

    

    X = SyntBlockColr(SytenyBlockInfo)
    Largest = 0 
    
    gdd = GenomeDiagram.Diagram('Diagram Name')
    for item in LocationData:
    
        #Find All Startings and endings
        AllStarts = [int(val[2]) for val in item[1:]] 
        AllEnds = [int(val[3]) for val in item[1:]] 
        StartSeq = (int(min(AllStarts)))
        EndSeq = (int(max(AllEnds)))



        #Initialize features 
        gdt_features = gdd.new_track(1, greytrack=False, name=[item[0][0]],
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
                    print(NewFeat, val)
                    Q = key
                    
            #print(NewFeat, Q)
            feature = SeqFeature(FeatureLocation(NewFeat[2], \
            NewFeat[3]), strand=None)
            gds_features.add_feature(feature, color=Q, name=str(NewFeat[1]), label=True)

    NameDraw = str(alignmentname) + "_genomediagram.pdf"
    gdd.draw(format='linear', pagesize=(15*cm,15*cm), fragments=1,
                         start=0, end=Largest)
    gdd.write(NameDraw, "pdf")
    
