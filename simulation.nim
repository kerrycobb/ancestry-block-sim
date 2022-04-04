import ./ancestry
import ./junctions
import ./haplotype
import ./deme
import std/random
# import std/io
import std/strutils

let
  generations = 1_000 
  demeSize = 1_000
  halfSize = int(demeSize/2)
  length = 1_000_000
  sitePos = 500_000
  recombRate = 1e-8
  recombProb = recombRate * length.float
  selection = rand(1.0)
  outfile = "simulation-output.txt"

var 
  startDeme = concat(
    newDemeWith(halfSize, singleSiteHaplotype(length, sitePos, ancSp1)),
    newDemeWith(halfSize, singleSiteHaplotype(length, sitePos, ancSp2)))
  endDeme = newDeme[Haplotype[SingleSite]](demeSize)

for gen in 0 ..< generations:
  endDeme = newDeme[Haplotype[SingleSite]](demeSize)
  while endDeme.size < endDeme.cap:
    let  
      a = sample(startDeme.members)
      b = sample(startDeme.members)
    if a.data.allele == b.data.allele:
      if recombProb < rand(1.0): 
        let bp = rand(length-1) # TODO: Decide if this should be length - 1 or just length 
        endDeme.members[endDeme.size] = recombine(a, b, bp) 
        endDeme.size.inc
      else:
        endDeme.members[endDeme.size] = a
        endDeme.size.inc
    else:
      if selection < rand(1.0): 
        if recombProb < rand(1.0): 
          let bp = rand(length-1) # TODO: Decide if this should be length - 1 or just length 
          endDeme.members[endDeme.size] = recombine(a, b, bp) 
          endDeme.size.inc
        else:
          endDeme.members[endDeme.size] = a
          endDeme.size.inc
  startDeme = endDeme 

var output = open(outfile, fmWrite)
for i in endDeme.members:
  # output.writeLine("")
  output.writeLine(i.junctions.join(","))
output.close()