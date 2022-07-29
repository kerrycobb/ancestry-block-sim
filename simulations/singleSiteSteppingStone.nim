## Simulate heterozygote disadvantage with a stepping stone model with demes of 
## randomly mating hermaphrodites with recombination

import ../src/ancestry
import ../src/junctions
import ../src/singleSiteHaplotype
import ../src/deme
import std/random
import std/strutils
import std/sequtils
import std/json

proc mating(startDeme: Deme, endDeme: var Deme, recombProb, selection: float): bool = 
  var
    hap1 = sample(startDeme.members) 
    hap2 = sample(startDeme.members) 
  if hap1.data.allele == hap2.data.allele: 
    if rand(1.0) < recombProb:
      endDeme.add(recombine(hap1, hap2, rand(hap1.len-2)))
    else:
      endDeme.add(hap1)
    result = true
  else:
    if rand(1.0) < selection:
      if rand(1.0) < recombProb:
        endDeme.add(recombine(hap1, hap2, rand(hap1.len-2)))
      else:
        endDeme.add(hap1)
      result = true
    else:
      result = false

proc migration(source: Deme, recipient: var Deme, nMigrants: int, recombProb, selection: float) = 
  var migCnt = 0
  while migCnt < nMigrants:
    var success = mating(source, recipient, recombProb, selection)
    if success:
      migCnt += 1

proc simulation(generations, demeNum, demeSize, length, position: int, recombRate: float, 
    seed: int, outfile: string) = 
  ## generations: number of generations
  ## demeNum: number of demes in stepping stone model
  ## demeSize: size of demes
  ## length: length of genomic region simulated
  ## position: position within genomic region under selection
  ## recombRate: per base probability of recombination occuring
  ## seed: starting seed
  ## outfile: output file path

  assert demeNum mod 2 == 0 # Deme number must be even
  assert position < length # Position must be smaller than the length
  randomize(seed)
  let
    recombProb = recombRate * (length).float
    selection = rand(1.0) # Probability that a heterozygote will be inherited
  var
    startDemes = newSeqWith(demeNum, newDeme[Haplotype[SingleSite]](demeSize))
    endDemes = newSeqWith(demeNum, newDeme[Haplotype[SingleSite]](demeSize))

  # Populate sp1 Demes
  for deme in 0 ..< demeNum div 2:
    startDemes[deme] = newDemeWith(demeSize, singleSiteHaplotype(length, position, ancSp1)) 
  
  # Populate sp2 Demes
  for deme in demeNum div 2 ..< demeNum:
    startDemes[deme] = newDemeWith(demeSize, singleSiteHaplotype(length, position, ancSp2)) 
  
  for gen in 0 ..< generations:
    # Clear endDeme to accept new haplotypes
    for i in 0 ..< demeNum:
      endDemes[i].empty()
  
    # Migration from source to edge
    endDemes[0].add(singleSiteHaplotype(length, position, ancSp1))
    endDemes[0].add(singleSiteHaplotype(length, position, ancSp1))
    endDemes[demeNum - 1].add(singleSiteHaplotype(length, position, ancSp2))
    endDemes[demeNum - 1].add(singleSiteHaplotype(length, position, ancSp2))
  
    # Migration from neighbor to edge
    migration(startDemes[1], endDemes[0], 2, recombProb, selection)
    migration(startDemes[0], endDemes[1], 2, recombProb, selection)
  
    # Migration from edge to neighbor
    migration(startDemes[demeNum-2], endDemes[demeNum-1], 2, recombProb, selection)
    migration(startDemes[demeNum-1], endDemes[demeNum-2], 2, recombProb, selection)
  
    # Migration between non edge demes
    for i in 1 ..< demeNum - 2:
      migration(startDemes[i], endDemes[i+1], 2, recombProb, selection)
      migration(startDemes[i+1], endDemes[i], 2, recombProb, selection)
  
    # Sample from within deme to populate next generation 
    for i in 0 ..< demeNum:
      while endDemes[i].size < endDemes[i].cap:
        discard mating(startDemes[i], endDemes[i], recombProb, selection)
    startDemes = endDemes
  
  var 
    output = open(outfile, fmWrite)
    metaData = %*{
      "seed": seed, 
      "selection": selection,
      "generations": generations, 
      "deme size": demeSize,
      "deme number": demeNum,
      "length": length,
      "site position": position,
      "recombination rate": recombRate,
    }
  output.writeLine($metaData)
  for i in endDemes:
    output.writeLine(i.members[0].junctions.join(","))
  output.close()

when isMainModule:
  import cligen; dispatch(simulation)