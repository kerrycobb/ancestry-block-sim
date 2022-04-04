import std/sequtils

type 
  Deme*[T] = object
    cap*: int
    size*: int
    members*: seq[T]

proc newDeme*[T](cap: int): Deme[T] = 
  Deme[T](cap:cap, size:0, members:newSeq[T](cap)) # TODO: Try newSeqOfCap for faster performance

template newDemeWith*(n: int, init: untyped): untyped = # naming the n paramter size to be conistent creates a bug 
  var result = Deme[typeof(init)](cap:n, size:n, members:newSeqWith(n, init))
  move(result)

proc concat*[T](demes: varargs[Deme[T]]): Deme[T] =
  for d in demes:
    result.cap += d.cap
    result.size += d.size
  result.members = concat(map(demes, proc(x: Deme[T]): seq[T] = x.members )) 

when isMainModule:
  var 
    d1 = newDemeWith(1, 1)
    d2 = newDemeWith(1, 2)
    c = concat(d1, d2)
  assert c.members == @[1, 2]

  echo "Tests Passing"

