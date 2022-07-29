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

proc mating(startDeme: Deme, endDeme: var Deme, recombRate, selection: float) = 
  var
    hap1 = sample(startDeme.members) 
    hap2 = sample(startDeme.members) 
  if hap1.data.allele == hap2.data.allele: 
    if rand(1.0) < recombRate:
      endDeme.add(recombine(hap1, hap2, rand(hap1.len-1)))
    else:
      endDeme.add(hap1)
  else:
    if rand(1.0) < selection:
      if rand(1.0) < recombRate:
        endDeme.add(recombine(hap1, hap2, rand(hap1.len-1)))
      else:
        endDeme.add(hap1)

proc simulation(replicates, generations, sampleSize, demeSize, length, 
    position: int, recombRate: float, seed: int, outfile: string) = 
  ## generations: number of generations
  ## sampleSize: number of samples to output from population 
  ## demeSize: size of deme
  ## length: length of genomic region simulated
  ## position: position within genomic region under selection
  ## recombRate: probability of a single recombination event per mating event 
  ## seed: starting seed
  ## outfile: output file path

  assert demeSize mod 2 == 0 # Deme size must be divisible by two 
  assert position < length # Position must be smaller than the length
  
  randomize(seed)
  let
    halfSize = demeSize div 2

  # Record some meta data in output file
  let 
    output = open(outfile, fmWrite)
    simulationMetaData = %*{
      "simulation type": "single site pancmictic",
      "replicates": replicates,
      "seed": seed, 
      "generations": generations, 
      "deme size": demeSize,
      "length": length,
      "site position": position,
      "recombination rate": recombRate,
    } 
  output.writeLine($simulationMetaData)

  # Simulation replicates
  for replicate in 0 ..< replicates:
    # record replicate meta data
    let 
      selection = rand(1.0)
      replicateMetaData = %*{
        "replicate": replicate,
        "selection": selection,
      }
    output.writeLine($replicateMetaData) 

    # Join demes with of two species
    var 
      startDeme = concat(
        newDemeWith(halfSize, singleSiteHaplotype(length, position, ancSp1)),
        newDemeWith(halfSize, singleSiteHaplotype(length, position, ancSp2)))
      endDeme = newDeme[Haplotype[SingleSite]](demeSize)

    # Simulate each generation 
    for gen in 0 ..< generations:
      endDeme.empty()
      while endDeme.size < endDeme.cap:
        mating(startDeme, endDeme, recombRate, selection)
      startDeme = endDeme 

    # Record replicate simulation data
    for i in 0..sampleSize-1:
      output.writeLine($i & ":" & endDeme.members[i].junctions.join(","))

  output.close()

when isMainModule:
  import cligen; dispatch(simulation)
