## This module contains the Junctions type wich is a `distinct seq[int]` with some
## borrowed procs. 
## Each int in Junctions corresponds to the beginning of a new ancestry block.
## Ancestry is assumed to be of species 1 until a junction is encountered.
## This means an emtpy Junctions object is assumed to have ancestry from only
## species 1. A junctions object with only a zero is assumed to have ancestry 
## only from species 2. Any subsequent int indicates a transition to the
## alternative state starting at that index. 
## Junctions of two haplotypes can be recombined with `recombine`. See doc string 
## for details.

import ./ancestry
import std/strutils

type Inheriting* = enum inheritingA, inheritingB

proc toggle*(i: var Inheriting) = 
  i = cast[Inheriting](i.ord xor inheritingB.ord)


type Junctions* = distinct seq[int] 

proc len*(j: Junctions): int {.borrow.}
proc `$`*(j: Junctions): string {.borrow.}
proc `==`*(a, b: Junctions): bool {.borrow.}
proc add*(j: var Junctions, n: int) {.borrow.}
proc `[]`*(a: Junctions, i: int): int = seq[int](a)[i]
proc join*(a: Junctions, s: string): string {.borrow.}

iterator items*(j: Junctions): int =  
  var i = 0
  while i < len(j):
    yield j[i]
    inc(i)

proc newJunctions*(sp: Ancestry): Junctions = 
  case sp
  of ancSp1:
    Junctions(@[])
  of ancSp2:
    Junctions(@[0])

proc recombine*(a, b: Junctions, breakPoints: seq[int]): Junctions = 
  ## Recombines the junctions of two haplotypes. The position of a breakPoint 
  ## corresponds to the right site of a recombination breakpoint and the beginning of 
  ## the new ancestry block being joined to. 
  ## Iters over existing junctions, adding them to new sequence of junctions if 
  ## they should be inherited. If breakpoint creates new junction it is added. 
  ## If breakpoint eliminates ancestry junction by bringing back together blocks
  ## that have the same ancestry, the junction is excluded.
  var 
    breakPoints = breakPoints
    aAncestryState = ancSp1
    bAncestryState = ancSp1
    inheriting = inheritingA 
    aix = 0
    bix = 0
  for bp in breakPoints:
    var
      aJunctionIsBreak = false 
      bJunctionIsBreak = false 
    # Iter over junctions up to breakpoint 
    case inheriting
    of inheritingA: 
      for i in aix ..< a.len: # Add inherited junctions
        if a[i] < bp:
          result.add(a[i])
          aAncestryState.toggle()
          aix.inc()
        else: 
          if a[i] == bp:
            aJunctionIsBreak = true
          break
      for i in bix ..< b.len: # Skip over not inherited junctions
        if b[i] < bp:
          bAncestryState.toggle()
          bix.inc()
        else: 
          if b[i] == bp:
            bJunctionIsBreak = true
          break
    of inheritingB: 
      for i in bix ..< b.len: # Add inherited junctions
        if b[bix] < bp:
          result.add(b[bix])
          bAncestryState.toggle()
          bix.inc()
        else:
          if b[i] == bp:
            bJunctionIsBreak = true
          break
      for i in aix ..< a.len: # Skip over non inherited junctions
        if a[i] < bp:
          aAncestryState.toggle()
          aix.inc()
        else: 
          if a[i] == bp:
            aJunctionIsBreak = true
          break
    # Determine if break point should be ancestry junction
    if aJunctionIsBreak and bJunctionIsBreak:
      if aAncestryState == bAncestryState: # Leave junction already shared by both haplotypes, will be skipped otherwise 
        result.add(bp)
      aix.inc()
      bix.inc()
      aAncestryState.toggle()
      bAncestryState.toggle()
    else: 
      if aAncestryState == bAncestryState:
        case inheriting 
        of inheritingA:
          if bJunctionIsBreak:
            result.add(bp)
            bix.inc()
            bAncestryState.toggle()
        of inheritingB: 
          if aJunctionIsBreak:
            result.add(bp)
            aix.inc()  
            aAncestryState.toggle()

      else: # aAncestryState != bAncestryState
        case inheriting
        of inheritingA: 
          if bJunctionIsBreak:
            bix.inc()
            bAncestryState.toggle()
          else:
            result.add(bp)
        of inheritingB: # If not inheritingFromA
          if aJunctionIsBreak:
            aix.inc()
            aAncestryState.toggle()
          else:
            result.add(bp)
    
    # Shift inheritance to other haplotype for next loop or remainder of junctions 
    inheriting.toggle()

  # Get remainder of breakpoints from inherited segment 
  case inheriting
  of inheritingA:
    for i in aix ..< a.len:
      result.add(a[i])
  of inheritingB:
    for i in bix ..< b.len:
      result.add(b[i])

proc recombine*(a, b: Junctions, breakPoint: int): Junctions = 
  recombine(a, b, @[breakPoint])

######################
# Original recombination proc that handles only a single breakpoint.
# Could resurrect this if speedup is necessary but it's less flexible and more 
# code that needs testing
#
# proc recombine(a, b: Haplotype, breakPoint: int): Haplotype = 
#   # Get breakpoints for first half of recombined haplotype for haplotype 1
#   var 
#     aAncestryState = ancestrySp1 
#     aJuncEqualBreak = false 
#   for i in 0 ..< a.len:
#     if a[i] < breakPoint: 
#       result.add(a[i])
#       aAncestryState.toggle
#     else:
#       if a[i] == breakPoint:
#         aJuncEqualBreak = true
#       break
#   # Skip over junctions up to new breakpoint for haplotype 2
#   var 
#     bAncestryState = ancestrySp1
#     bix = 0
#     bJuncEqualBreak = false
#   for i in 0 ..< b.len:
#     if b[i] < breakPoint:
#       bix.inc
#       bAncestryState.toggle
#     else:
#       if b[i] == breakPoint:
#         bJuncEqualBreak = true 
#       break
#   # Determine if break point should be ancestry junction
#   if aJuncEqualBreak and bJuncEqualBreak:
#     if aAncestryState == bAncestryState: 
#       result.add(breakPoint)
#     bix.inc 
#   else:
#     if aAncestryState != bAncestryState:
#       if not bJuncEqualBreak:
#         result.add(breakPoint)
#       else:
#         bix.inc 
#     else:
#       if bJuncEqualBreak:
#         result.add(breakPoint)
#         bix.inc 
#   # Get remainder of breakpoints from haplotype b 
#   for i in bix ..< b.len:
#     result.add(b[i])

when isMainModule:
  # Test ancestry toggle
  var a = ancSp1
  a.toggle
  assert a == ancSp2 
  a.toggle
  assert a == ancSp1

  # Test newJunctions  
  assert newJunctions(ancSp1) == Junctions(@[]) 
  assert newJunctions(ancSp2) == Junctions(@[0]) 

  ## Tests with single breakpoint
  # Breakpoint at zero
  assert recombine(Junctions(@[]),  Junctions(@[ ]), 0) == Junctions(@[ ])
  assert recombine(Junctions(@[0]), Junctions(@[0]), 0) == Junctions(@[0])
  assert recombine(Junctions(@[ ]), Junctions(@[0]), 0) == Junctions(@[0])
  assert recombine(Junctions(@[0]), Junctions(@[ ]), 0) == Junctions(@[ ])
  # Tests with no existing junctions
  assert recombine(Junctions(@[0]), Junctions(@[0]), 50) == Junctions(@[0    ])
  assert recombine(Junctions(@[ ]), Junctions(@[ ]), 50) == Junctions(@[     ])
  assert recombine(Junctions(@[0]), Junctions(@[ ]), 50) == Junctions(@[0, 50])
  assert recombine(Junctions(@[ ]), Junctions(@[0]), 50) == Junctions(@[50   ])
  # Test when break point corresponds with junction
  assert recombine(Junctions(@[0, 50]), Junctions(@[     ]), 50) == Junctions(@[0, 50])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0    ]), 50) == Junctions(@[0    ])
  assert recombine(Junctions(@[0, 50]), Junctions(@[50   ]), 50) == Junctions(@[0    ])
  assert recombine(Junctions(@[     ]), Junctions(@[0, 50]), 50) == Junctions(@[     ])
  assert recombine(Junctions(@[0    ]), Junctions(@[0, 50]), 50) == Junctions(@[0, 50])
  assert recombine(Junctions(@[50   ]), Junctions(@[0, 50]), 50) == Junctions(@[     ])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0, 50]), 50) == Junctions(@[0, 50])
  # Same as above but with junctions after breakpoint 
  assert recombine(Junctions(@[0, 50]), Junctions(@[       60]), 50) == Junctions(@[0, 50, 60])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0    , 60]), 50) == Junctions(@[0    , 60])
  assert recombine(Junctions(@[0, 50]), Junctions(@[50   , 60]), 50) == Junctions(@[0    , 60])
  assert recombine(Junctions(@[     ]), Junctions(@[0, 50, 60]), 50) == Junctions(@[       60])
  assert recombine(Junctions(@[0    ]), Junctions(@[0, 50, 60]), 50) == Junctions(@[0, 50, 60])
  assert recombine(Junctions(@[50   ]), Junctions(@[0, 50, 60]), 50) == Junctions(@[       60])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0, 50, 60]), 50) == Junctions(@[0, 50, 60])
  # Test with existing junctions not overlapping breakpoints 
  assert recombine(Junctions(@[0, 50]), Junctions(@[0, 50]), 25) == Junctions(@[0,  50])
  assert recombine(Junctions(@[50   ]), Junctions(@[0, 50]), 25) == Junctions(@[25, 50])
  assert recombine(Junctions(@[0    ]), Junctions(@[0, 50]), 25) == Junctions(@[0,  50])
  assert recombine(Junctions(@[     ]), Junctions(@[0, 50]), 25) == Junctions(@[25, 50])
  ## Tests with multiple breakpoints
  # With no existing junctions
  assert recombine(Junctions(@[0]), Junctions(@[0]), @[50, 100]) == Junctions(@[0    ])
  assert recombine(Junctions(@[ ]), Junctions(@[ ]), @[50, 100]) == Junctions(@[     ])
  assert recombine(Junctions(@[0]), Junctions(@[ ]), @[50, 100]) == Junctions(@[0, 50, 100])
  assert recombine(Junctions(@[ ]), Junctions(@[0]), @[50, 100]) == Junctions(@[50, 100   ])
  # With junctions that overlap breakpoint
  assert recombine(Junctions(@[0, 50]), Junctions(@[     ]), @[50, 100]) == Junctions(@[0, 50     ])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0    ]), @[50, 100]) == Junctions(@[0, 100    ])
  assert recombine(Junctions(@[0, 50]), Junctions(@[50   ]), @[50, 100]) == Junctions(@[0, 100    ])
  assert recombine(Junctions(@[     ]), Junctions(@[0, 50]), @[50, 100]) == Junctions(@[          ])
  assert recombine(Junctions(@[0    ]), Junctions(@[0, 50]), @[50, 100]) == Junctions(@[0, 50, 100])
  assert recombine(Junctions(@[50   ]), Junctions(@[0, 50]), @[50, 100]) == Junctions(@[100       ])
  assert recombine(Junctions(@[0, 50]), Junctions(@[0, 50]), @[50, 100]) == Junctions(@[0, 50     ])
  # TODO: Add more tests with multiple recombination events

  echo "junctions.nim Tests Passing"