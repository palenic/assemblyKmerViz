"""Set of functions to evaluate k-mer multiplicity in each position of a genome assembly
starting from meryl-lookup results.
"""

import pandas as pd
import random

def merylLookupWigToBed(input_wig,output_bed,verbose=True):
    """ Summarizes variableStep .wig files produced by 
    meryl-lookup -wig-count showing the multiplicity of the kmer
    starting at each position in the genome assembly.
    Produces a tab-separated .bed file that aggregates adjacent positions with 
    the same multiplicity. Scores >4 are aggregated as 5.
    The file is 0-based, half-open [start-1,end) and has fields:
    chrom	start	end .	multiplicity    .   start   start   color
    
    Parameters
    ----------
    input_wig: string 
        path to input wig file
    output_bed: string 
        path to output file
    verbose: boolean 
        print verbose output
    
    Author: Chiara Paleni
    """
    #colors from merqury
    gray = "0,0,0"
    red = "228,26,28"
    blue = "55,126,184" # light blue = "#56B4E9"
    green = "77,175,74"
    purple = "152,78,163"  # purple = "#CC79A7"
    orange = "255,127,0"  # orange = "#E69F00"
    yellow = "255,255,51"

    merqury_col = [gray, red, blue, green, purple, orange]
    
    o=open(output_bed,"w")
    
    with input_wig as f:
    
        # read first line to check it has the correct headers
        # correct header looks like this:
        # variableStep chrom=tigName
        
        wigLine=f.readline().strip().split("\t")
        if(not wigLine[0].startswith("variableStep")):
            raise TypeError("Missing variableStep header. Check your input file.")
        tigname=wigLine[0].split(" ")[1].split("=")[1]
        if(verbose):
            print("contig",tigname)
        
        startPos=0
        endPos=0
        score=-1
        strand="."
        name="."
        
        
        #score is the multiplicity. initialized as -1 because it cannot have this value in real lines
        #start reading following lines
        
        for x in f:
            wigLine=x.strip().split("\t")
            #print(wigLine)
            if(wigLine[0].startswith("variableStep")):
            
                # if i am at the start of a new chromosome, write previous region to file 
                # and initialize tigName, startPos, endPos, score
                
                newBedLine=[tigname,startPos,endPos,name,score,strand,startPos,startPos,merqury_col[score]]
                if(score!=-1):
                    o.write("\t".join(str(item) for item in newBedLine))
                    o.write("\n")
                
                tigname=wigLine[0].split(" ")[1].split("=")[1]
                if(verbose):
                    print("contig",tigname)
                
                startPos=0
                endPos=0
                score=-1
                
            else:
                #if i am in a "normal" line, two cases:                
                #the score is different: write previous region to file and create a new one.
                #the score is the same as the previous region: aggregate.
                #print(int(a[1]))
                newScore=int(wigLine[1])
                if(newScore>=5):
                    newScore=5
                    #i will keep track for 1,2,3,4,>4 like merqury does.
                
                if(score!=int(newScore)):
                    #print(score)
                    newBedLine=[tigname,startPos,endPos,name,score,strand,startPos,startPos,merqury_col[score]]
                    
                    if(score!=-1):
                        # score is -1 if i have just started reading a new chromosome.
                        # so, even if the newScore is different from the old one,
                        # i don't need to write the region in that case since it is empty.
                        #print(newBedLine)
                        o.write("\t".join(str(item) for item in newBedLine))
                        o.write("\n")
                    
                    # the new region will have the following coordinates.    
                    endPos+=1
                    startPos=int(wigLine[0])-1
                    score=newScore                    
                
                else:
                    # extend previous region.
                    endPos+=1
                    #print(wigLine)
    #for the last line:
    newBedLine=[tigname,startPos,endPos,name,score,strand,startPos,startPos,merqury_col[score]]
    if(score!=-1):
        #print(newBedLine)
        o.write("\t".join(str(item) for item in newBedLine))
        o.write("\n")
    o.close()
    if(verbose):
        print("done!")

def merylLookupBedProcessing(input_bed,output_csv,relativeCounts=False,verbose=True):
    """ Takes a bed file produced by merylLookupWigToBed
    and turns it into a .csv file containing for
    each contig its name, length, and counts of how many
    positions in the contig have that multiplicity,
    optionally scaled based on contig length.
    Scores >4 are aggregated as 5.
    Produces a .csv file with columns: contigName,contigLength,
    counts_0,counts_1,counts_2,counts_3,counts_4,counts_5_or_more.
    The file is sorted by contigName.
    
    Parameters
    ----------
    input_bed: string 
        path to input bed file
    output_csv: string 
        path to output csv file
    relativeCounts: boolean
        if true, returns the percentage of bases with that
        multiplicity instead of the counts (default False)
    verbose: boolean 
        print verbose output
        
    Author: Chiara Paleni
    """
    counts=[]
    with open(input_bed) as f:
        tigCounts=[0,0,0,0,0,0]
        tigName=None
        tigNames=[]
        for x in f:
            l=x.strip().split()
            #print(l)
            if(tigName!=l[0]):
                if(tigName is not None):
                    #print(tigCounts)
                    counts.append(tigCounts)
                tigName=l[0]
                tigCounts=[0,0,0,0,0,0]
                tigNames.append(tigName)
                if(verbose):
                    print("found contig ",tigName)
            mult=int(l[4])
            #print(mult)
            length=int(l[2])-int(l[1])
            tigCounts[mult]+=length
        #print(tigCounts)
        counts.append(tigCounts)
    #print(counts)
    tiglen=[sum(x) for x in counts]
    counts_df=pd.DataFrame(counts,columns=[0,1,2,3,4,5])
    counts_df["contig"]=pd.Series(tigNames)
    counts_df["tigLen"]=pd.Series(tiglen)
    #counts_df.columns
    counts_df=counts_df.reindex(columns=['contig','tigLen',0, 1, 2, 3, 4, 5])
    counts_df = counts_df.sort_values(['contig'], ascending=True)
    counts_df=counts_df.reindex()
    counts_df=counts_df.reset_index(drop=True)
    if(relativeCounts):
        for i in range(6):
            counts_df[i]=counts_df[i]/counts_df["tigLen"]
    if(verbose):
        print("done!")
    counts_df.to_csv(output_csv,index=False)
    
    
    
    
def merylLookupWigToBedBinned(input_wig,output_bed,bin_width=100,verbose=True):
    """ Summarizes variableStep .wig files produced by 
    meryl-lookup -wig-count showing the multiplicity of the kmer
    starting at each position in the genome assembly.
    Produces a tab-separated .bed file in which multiplicity inside bins 
    of the specified width (default -b 100) is computed with majority vote 
    of all positions inside the bin. Scores >4 are aggregated as 5.
    The file is 0-based, half-open [start-1,end) and has fields:
    chrom	start	end .	multiplicity    .   start   start   color
    
    Parameters
    ----------
    input_wig: string 
        path to input wig file
    output_bed: string 
        path to output file
    bin_width: int (default 2)
        width of bins for reporting the multiplicity. The score reported for 
        each region is decided through majority vote of all positions 
        belonging to that interval.
    verbose: boolean (default True)
        print verbose output
    
    Author: Chiara Paleni
    """
    #colors from merqury
    gray = "0,0,0"
    red = "228,26,28"
    blue = "55,126,184" # light blue = "#56B4E9"
    green = "77,175,74"
    purple = "152,78,163"  # purple = "#CC79A7"
    orange = "255,127,0"  # orange = "#E69F00"
    yellow = "255,255,51"

    merqury_col = [gray, red, blue, green, purple, orange]
    
    def majorityVote(votes):
        maxScore=0
        maxPos=[]
        for i in range(0,len(votes)):
            if(votes[i]>maxScore):
                maxScore=votes[i]
                maxPos=[]
                maxPos.append(i)
            elif (votes[i]==maxScore):
                maxPos.append(i)
        if(len(maxPos)>1):
            maxPos=random.choice(maxPos)
        else:
            maxPos=maxPos[0]
        return maxPos
    o=open(output_bed,"w")
    
    #def writeBedLine
    
    with input_wig as f:
    
        # read first line to check it has the correct headers
        # correct header looks like this:
        # variableStep chrom=tigName
        
        wigLine=f.readline().strip().split("\t")
        if(not wigLine[0].startswith("variableStep")):
            raise TypeError("Missing variableStep header. Check your input file.")
        tigname=wigLine[0].split(" ")[1].split("=")[1]
        if(verbose):
            print("contig",tigname)
        
        startPos=0
        endPos=startPos+bin_width
        newPos=1
        #score=-1
        votes=[0,0,0,0,0,0]
        strand="."
        name="."
        
        
        #score is the multiplicity. initialized as -1 because it cannot have this value in real lines
        #start reading following lines
        
        for x in f:
            wigLine=x.strip().split("\t")
            #print(wigLine)
            if(wigLine[0].startswith("variableStep")):
            
                # if i am at the start of a new chromosome, write previous region to file 
                # and initialize tigName, startPos, endPos, score
                # if i am at the start of the file (score=-1) do nothing
                
                if(votes!=[0,0,0,0,0,0]):
                    score=majorityVote(votes)
                    newBedLine=[tigname,startPos,newPos,name,score,strand,startPos,startPos,merqury_col[score]]
                    o.write("\t".join(str(item) for item in newBedLine))
                    o.write("\n")
                    #print(newBedLine)
                
                tigname=wigLine[0].split(" ")[1].split("=")[1]
                if(verbose):
                    print("contig",tigname)
                
                startPos=0
                endPos=startPos+bin_width
                votes=[0,0,0,0,0,0]
                
            else:
                #if i am in a "normal" line, two cases:                
                #i am inside the bin: add score to votes
                #i am outside the bin: vote and write previous bin, then 
                #    initialize new one
                newPos=int(wigLine[0])
                newScore=int(wigLine[1])
                score=newScore
                if(newScore>=5):
                    newScore=5
                    #i will keep track for 1,2,3,4,>4 like merqury does.
                if(newPos>endPos):
                    score=majorityVote(votes)
                    newBedLine=[tigname,startPos,endPos,name,score,strand,startPos,startPos,merqury_col[score]]
                    o.write("\t".join(str(item) for item in newBedLine))
                    o.write("\n")
                    startPos=endPos
                    endPos=startPos+bin_width
                    votes=[0,0,0,0,0,0]
                    votes[newScore]+=1
                    #print(newBedLine)
                    #print(votes)
                else:
                    votes[newScore]+=1
                    #print(votes)
    #for the last line:
    
    if(votes!=[0,0,0,0,0,0]):
        score=majorityVote(votes)
        newBedLine=[tigname,startPos,newPos,name,score,strand,startPos,startPos,merqury_col[score]]
        #print(newBedLine)
        o.write("\t".join(str(item) for item in newBedLine))
        o.write("\n")
    o.close()
    if(verbose):
        print("done!")
