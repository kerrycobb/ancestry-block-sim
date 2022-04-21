
import ./ancestry
import ./junctions
import std/sequtils

type 
  Allele* = enum allele1, allele2

  SingleSite* = object
    pos*: int
    allele*: Allele

  Haplotype*[T] = object
    len*: int 
    junctions*: Junctions
    data*: T

proc toggle*(a: var Allele) = 
  a =  cast[Allele](a.ord xor allele2.ord)

proc singleSiteHaplotype*(len, pos: int, sp: Ancestry): Haplotype[SingleSite] =
  result.len = len
  result.junctions = newJunctions(sp)
  case sp 
  of ancSp1:
    result.data = SingleSite(pos:pos, allele:allele1)
  of ancSp2:
    result.data = SingleSite(pos:pos, allele:allele2)

proc recombine*(a, b: Haplotype[SingleSite], breakPoints: seq[int]): Haplotype[SingleSite] =
  assert a.len == b.len
  assert a.data.pos == a.data.pos
  apply(breakPoints, proc(x: int) = assert x < a.len)
  apply(breakPoints, proc(x: int) = assert x > 0)
  result.len = a.len
  result.data.pos = a.data.pos
  result.junctions = recombine(a.junctions, b.junctions, breakPoints)
  for junc in result.junctions:
    if junc <= a.data.pos: 
      result.data.allele.toggle
    else:
      break

proc recombine*(a, b: Haplotype[SingleSite], breakPoint: int): Haplotype[SingleSite] =
  recombine(a, b, @[breakPoint])


when isMainModule:
  var a = allele1
  a.toggle
  assert a == allele2 
  a.toggle
  assert a == allele1

  assert singleSiteHaplotype(5, 5, ancSp1) == HaploType[SingleSite](len:5, junctions: Junctions(@[ ]), data: SingleSite(pos:5, allele:allele1))
  assert singleSiteHaplotype(5, 5, ancSp2) == HaploType[SingleSite](len:5, junctions: Junctions(@[0]), data: SingleSite(pos:5, allele:allele2))

  assert recombine(singleSiteHaplotype(5, 5, ancSp1), singleSiteHaplotype(5, 5, ancSp2), 3) == Haplotype[SingleSite](len:5, junctions: Junctions(@[3   ]), data: SingleSite(pos:5, allele:allele2))  
  assert recombine(singleSiteHaplotype(5, 5, ancSp2), singleSiteHaplotype(5, 5, ancSp1), 3) == Haplotype[SingleSite](len:5, junctions: Junctions(@[0, 3]), data: SingleSite(pos:5, allele: allele1))  
  
  echo "singleSiteHaplotype.nim Tests Passing"



# TODO: Experiment with static int types for haplotype and allele, only useful if 
# numbers can be specified at compile time. Probably not something worth pursuing.

# import ./ancestry
# import ./junctions

# type 
#   Allele* = enum allele1, allele2

#   SingleSite*[Position: static[int]] = object
#     allele*: Allele

#   Haplotype*[Length: static[int], T] = object
#     junctions*: Junctions
#     data*: T

# proc toggle*(a: var Allele) = 
#   a =  cast[Allele](a.ord xor allele2.ord)

# proc singleSiteHaplotype*(len, pos: static[int], sp: Ancestry): Haplotype[len, SingleSite[pos]] =
#   result.junctions = newJunctions(sp)
#   case sp 
#   of ancSp1:
#     result.data = SingleSite[pos](allele: allele1)
#   of ancSp2:
#     result.data = SingleSite[pos](allele: allele2)

# proc recombine*[Length, Position](a, b: Haplotype[Length, SingleSite[Position]], 
#     breakPoints: seq[int]): Haplotype[Length, SingleSite[Position]] =
#   for i in breakPoints: doAssert i <= Length
#   result.junctions = recombine(a.junctions, b.junctions, breakPoints)
#   for junc in result.junctions:
#     if junc <= Position: 
#       result.data.allele.toggle
#     else:
#       break

# proc recombine*[Length, Position](a, b: Haplotype[Length, SingleSite[Position]], 
#     breakPoint: int): Haplotype[Length, SingleSite[Position]] =
#   recombine(a, b, @[breakPoint])

# if isMainModule:
#   var a = allele1
#   a.toggle
#   assert a == allele2 
#   a.toggle
#   assert a == allele1

#   assert singleSiteHaplotype(5, 5, ancSp1) == HaploType[5, SingleSite[5]](junctions: Junctions(@[]), data: SingleSite[5](allele: allele1))
#   assert singleSiteHaplotype(5, 5, ancSp2) == HaploType[5, SingleSite[5]](junctions: Junctions(@[0]), data: SingleSite[5](allele: allele2))

#   assert recombine(singleSiteHaplotype(5, 5, ancSp1), singleSiteHaplotype(5, 5, ancSp2), 3) == Haplotype[5, SingleSite[5]](junctions: Junctions(@[3   ]), data: SingleSite[5](allele: allele2))  
#   assert recombine(singleSiteHaplotype(5, 5, ancSp2), singleSiteHaplotype(5, 5, ancSp1), 3) == Haplotype[5, SingleSite[5]](junctions: Junctions(@[0, 3]), data: SingleSite[5](allele: allele1))  
  
#   echo "Tests Pass"



