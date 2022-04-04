## Ancestry type used to indicate the ancestral species of origin for a haplotype

type
  Ancestry* = enum ancSp1, ancSp2

proc toggle*(a: var Ancestry) = 
  ## Switch ancestry to opposite state
  a =  cast[Ancestry](a.ord xor ancSp2.ord)

when isMainModule:
  var a = ancSp1
  a.toggle
  assert a == ancSp2
  a.toggle
  assert a == ancSp1
  echo "Tests Passing"