import ../src/ancestry
import ../src/junctions
import ../src/singleSiteHaplotype
import ../src/deme
import std/random
import std/strutils
import std/sequtils

let
  generations = 100 
  demeSize = 100
  halfSize = demeSize div 2
  length = 1_000_000
  sitePos = 500_000
  recombRate = 1e-8
  recombProb = recombRate * (length-1).float
  selection = rand(1.0)
  outfile = "singleSitePanmictic-output.txt"

var 
  startDeme = concat(
    newDemeWith(halfSize, singleSiteHaplotype(length, sitePos, ancSp1)),
    newDemeWith(halfSize, singleSiteHaplotype(length, sitePos, ancSp2)))
  endDeme = newDeme[Haplotype[SingleSite]](demeSize)

for gen in 0 ..< generations:
  endDeme.empty()
  while endDeme.size < endDeme.cap:
    let  
      hap1 = sample(startDeme.members)
      hap2 = sample(startDeme.members)
    if hap1.data.allele == hap2.data.allele:
      if recombProb < rand(1.0): 
        endDeme.add(recombine(hap1, hap2, rand(1..hap1.len-1))) 
      else:
        endDeme.add(hap1)
    else:
      if selection < rand(1.0): 
        if recombProb < rand(1.0): 
          let bp = rand(length-1)
          endDeme.add(recombine(hap1, hap2, rand(1..hap1.len-1)))
        else:
          endDeme.add(hap1)
  startDeme = endDeme 

var output = open(outfile, fmWrite)
for i in endDeme.members:
  output.writeLine(i.junctions.join(","))
output.close()


