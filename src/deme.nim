import std/sequtils

type 
  Deme*[T] = object
    cap*: int
    size*: int
    members*: seq[T]

proc newDeme*[T](cap: int): Deme[T] = 
  Deme[T](cap:cap, size:0, members:newSeqOfCap[T](cap)) # TODO: Try newSeqOfCap for faster performance

template newDemeWith*(n: int, init: untyped): untyped = # naming the n paramter size to be conistent creates a bug 
  var result = Deme[typeof(init)](cap:n, size:n, members:newSeqWith(n, init))
  move(result)

proc `[]`*[T](deme: Deme[T], i: int): T = 
  deme.members[i]

proc `[]=`*[T](deme: var Deme[T], i: int, v: T) = 
  deme.members[i] = v

proc concat*[T](demes: varargs[Deme[T]]): Deme[T] =
  for d in demes:
    result.cap += d.cap
    result.size += d.size
  result.members = concat(map(demes, proc(x: Deme[T]): seq[T] = x.members )) 

proc add*[T](deme: var Deme[T], member: T) = 
  assert deme.cap != deme.size
  deme.members.setLen(deme.size + 1)
  deme.members[deme.size] = member 
  deme.size += 1

proc isFull*[T](deme: Deme[T]): bool = 
  if deme.size == deme.cap:
    result = true

proc empty*[T](deme: var Deme[T]) = 
  deme.members.setLen(0)
  deme.size = 0

when isMainModule:
  var 
    d1 = newDemeWith(1, 1)
    d2 = newDemeWith(1, 2)
    c = concat(d1, d2)
  assert c.members == @[1, 2]
  c.empty()
  c.add(1)
  assert c.members == @[1]

  echo "deme.nim Tests Passing"