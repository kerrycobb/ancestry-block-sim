## Simulate heterozygote disadvantage in a single panmictic population of 
## randomly mating hermaphrodites with recombination

import ../src/ancestry
import ../src/junctions
import ../src/singleSiteHaplotype
import ../src/deme
import std/random
import std/strutils
import std/sequtils
import std/json

proc mating(startDeme: Deme, endDeme: var Deme, recombProb, selection: float) = 
  var
    hap1 = sample(startDeme.members) 
    hap2 = sample(startDeme.members) 
  if hap1.data.allele == hap2.data.allele: 
    if rand(1.0) < recombProb:
      endDeme.add(recombine(hap1, hap2, rand(hap1.len-2)))
    else:
      endDeme.add(hap1)
  else:
    if rand(1.0) < selection:
      if rand(1.0) < recombProb:
        endDeme.add(recombine(hap1, hap2, rand(hap1.len-2)))
      else:
        endDeme.add(hap1)

proc simulation(generations, sampleSize, demeSize, length, position: int, recombRate: float, 
    seed: int, outfile: string) = 
  ## generations: number of generations
  ## demeSize: size of deme
  ## length: length of genomic region simulated
  ## position: position within genomic region under selection
  ## recombRate: per base probability of recombination occuring
  ## seed: starting seed
  ## outfile: output file path

  assert demeSize mod 2 == 0 # Deme size must be even
  assert position < length # Position must be smaller than the length
  
  randomize(seed)
  let
    halfSize = demeSize div 2
    recombProb = recombRate * (length).float
    selection = rand(1.0)
  
  var 
    startDeme = concat(
      newDemeWith(halfSize, singleSiteHaplotype(length, position, ancSp1)),
      newDemeWith(halfSize, singleSiteHaplotype(length, position, ancSp2)))
    endDeme = newDeme[Haplotype[SingleSite]](demeSize)
  
  for gen in 0 ..< generations:
    endDeme.empty()
    while endDeme.size < endDeme.cap:
      mating(startDeme, endDeme, recombProb, selection)
    startDeme = endDeme 

  var 
    output = open(outfile, fmWrite)
    metaData = %*{
      "seed": seed, 
      "selection": selection,
      "generations": generations, 
      "deme size": demeSize,
      "length": length,
      "site position": position,
      "recombination rate": recombRate,
    } 
  output.writeLine($metaData)
  for i in 0..sampleSize-1:
    output.writeLine(endDeme.members[i].junctions.join(","))
  output.close()

when isMainModule:
  import cligen; dispatch(simulation)
