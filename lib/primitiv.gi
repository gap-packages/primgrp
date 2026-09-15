#############################################################################
##
#W  primitiv.gi       GAP primitive groups library          Alexander Hulpke
##
##
#Y  Copyright (C)  1999, School Math.&Comp. Sci., University of St Andrews
##
##  This file contains the routines for the primitive groups library
##

Unbind(PRIMGRP);

#############################################################################
##
#V  PRIMGRP
##  Generators, names and properties of the primitive groups.
##  entries are
##  1: id
##  2: size
##  3: Simple+2*Solvable
##  4: ONan-Scott-type
##  5: Collected suborbits
##  6: Transitivity
##  7: name
##  8: socle type
##  9: generators
BindGlobal("PRIMGRP", []);

#############################################################################
##
#F  PGAlt( <deg>, <nr> ) . . . . . . . . . . . the natural A_n and S_n entries
#F  PGSym( <deg>, <nr> )
##
##  Everything an entry records about the natural alternating and symmetric
##  group follows from the degree: the order, that A_n is simple and S_n is
##  not, the O'Nan-Scott type, the single suborbit of a 2-transitive group,
##  the transitivity, the name and the socle.  So the entry can be the call
##  itself.
##
##  A data file names these rather than calling them, as ["Alt"] and ["Sym"];
##  PRIMGrp evaluates that on first use.  Factorial(deg) is then not computed
##  for every degree in the file merely because the file was read.
##
BindGlobal("PGAlt",function(deg,nr)
  return [ nr, Factorial(deg)/2, 1, "2", [[deg-1,1]], deg-2,
           Concatenation("A(",String(deg),")"), ["A",deg,1], "Alt" ];
end);

BindGlobal("PGSym",function(deg,nr)
  return [ nr, Factorial(deg), 0, "2", [[deg-1,1]], deg,
           Concatenation("S(",String(deg),")"), ["A",deg,1], "Sym" ];
end);

#############################################################################
##
#F  PGPrime( <d> ) . . . . . . a subgroup of AGL(1,p) of order p*<d>, p prime
##
##  For prime degree p the affine primitive groups are exactly the subgroups of
##  AGL(1,p) that contain the translations, one for each divisor d of p-1.
##  Everything the entry records follows from p and d, so a data file stores
##  only d, as ["Prime",d].
##
##  The multiplicative group of GF(p) is cyclic, so its subgroup of order d is
##  unique: any element of order d generates it, and the group does not depend
##  on which primitive root PrimitiveRootMod happens to return.
##
##  The names are the ones the library already uses at these degrees: C(p) for
##  the translations alone, D(2*p) when d = 2, AGL(1, p) when d = p-1, and p:d
##  otherwise.
##
##  This starts at p = 5.  AGL(1,2) and AGL(1,3) are Sym(2) and Sym(3), which
##  the library names S(2), A(3) and S(3) and PGAlt and PGSym describe; and
##  Sym(3) is 3-transitive where AGL(1,p) is 2-transitive, so the rules below
##  would get those entries wrong rather than merely name them differently.
##
BindGlobal("PGPrime",function(d)
  return function(deg,nr)
    local a,name,gens,flags,trans;
    if deg < 5 then
      Error("PGPrime: AGL(1,",deg,") is Sym(",deg,
            "), which PGAlt and PGSym describe");
    fi;
    if (deg-1) mod d <> 0 then
      Error("PGPrime(",d,") at degree ",deg,": ",d," does not divide ",deg-1);
    fi;
    if d = 1 then
      name:=Concatenation("C(",String(deg),")");
      gens:=[];
    else
      a:=PowerModInt(PrimitiveRootMod(deg),(deg-1)/d,deg);
      gens:=[ [ [ a*Z(deg)^0 ] ] ];
      if d = 2 then
        name:=Concatenation("D(2*",String(deg),")");
      elif d = deg-1 then
        name:=Concatenation("AGL(1, ",String(deg),")");
      else
        name:=Concatenation(String(deg),":",String(d));
      fi;
    fi;
    if d = 1 then flags:=3; else flags:=2; fi;      # simple and solvable
    if d = deg-1 then trans:=2; else trans:=1; fi;  # AGL(1,p) is 2-transitive
    return [ nr, deg*d, flags, "1", [[d,(deg-1)/d]],
             trans, name, ["Z",deg,1], gens ];
  end;
end);

#############################################################################
##
#F  PRIMGRP_JohnsonSuborbits( <n>, <k> ) . . . . Alt(n) and Sym(n) on k-sets
##
##  The subdegrees of the Johnson scheme.  The stabiliser of a k-set S has one
##  orbit for each i = 1..k, holding the k-sets that agree with S in all but i
##  of its points: Binomial(k,i) choices of which points of S to drop, and
##  Binomial(n-k,i) of what to put in their place.  At n = 45, k = 2 that is
##  2*43 = 86 sets sharing a point with S, and 1*903 = 903 disjoint from it.
##
BindGlobal("PRIMGRP_JohnsonSuborbits",function(n,k)
  return Set(Collected(List([1..k],i->Binomial(k,i)*Binomial(n-k,i))));
end);

#############################################################################
##
#F  PGOnPointsGroup( <inner> ) . . . . PSL extended by words, on the points
##
##  <inner> is ["PSL",dim,q,<words>,<name>]: PSL(dim,q) on the points of
##  PG(dim-1,q), extended by one element for each word [i,j] in <words>,
##  diag(Z(q)^i,1,...,1) followed by the j-th power of the Frobenius.  That
##  reaches every group between PSL and PGammaL, including those the four named
##  ones miss, such as PSL(2,25).2_3 = PSL(2,25).<delta phi>.  The points are
##  the sorted normed vectors.
##
##  Modulo PSL, delta = diag(Z(q),1,...,1) has order g = gcd(dim,q-1) and the
##  Frobenius phi order f, where q = p^f.  First delta^i1 phi^j1 then
##  delta^i2 phi^j2 is delta^(i1 + i2*p^(f-j1)) phi^(j1+j2), and
##  PRIMGRP_ProjectiveExtension closes the words under that.
##
BindGlobal("PRIMGRP_ProjectiveExtension",function(dim,q,words)
  local p,f,g,mul,H,new,x,y,z;
  p:=Characteristic(GF(q));
  f:=Length(Factors(q));
  g:=Gcd(dim,q-1);
  mul:=function(a,b)
    return [ (a[1]+b[1]*PowerModInt(p,(f-a[2]) mod f,g)) mod g,
             (a[2]+b[2]) mod f ];
  end;
  H:=[[0,0]];
  new:=[[0,0]];
  while new <> [] do
    x:=Remove(new);
    for y in words do
      z:=mul(x,[y[1] mod g,y[2] mod f]);
      if not z in H then
        Add(H,z);
        Add(new,z);
      fi;
    od;
  od;
  return Set(H);
end);

BindGlobal("PGOnPointsGroup",function(inner)
  local dim,q,p,vecs,gens,w;
  dim:=inner[2];
  q:=inner[3];
  p:=Characteristic(GF(q));
  vecs:=Set(NormedRowVectors(GF(q)^dim));
  gens:=List(GeneratorsOfGroup(SL(dim,q)),m->Permutation(m,vecs,OnLines));
  for w in inner[4] do
    Add(gens,Permutation(
      DiagonalMat(Concatenation([Z(q)^w[1]],List([2..dim],i->One(GF(q))))),
      vecs,function(v,m) return List(OnLines(v,m),c->c^(p^w[2])); end));
  od;
  return GroupWithGenerators(gens);
end);

#############################################################################
##
#F  PGOnSetsGroup( <inner>, <k> ) . . . . . . an inner group on the k-subsets
##
##  Build the group an entry's field 9 asks for: PRIMGRP_InnerGroup turns the
##  description <inner> into a group, and PGOnSetsGroup returns its action on
##  the k-element subsets of its points.
##
##  The points are all k-subsets of [1..n] in sorted order, and the image of a
##  k-set is located by its rank in that order rather than searched for.  That
##  order is fixed; one given by Orbit would not be, and which permutation group
##  comes out depends on it.  The result is transitive only if the group is
##  k-homogeneous, so that is checked.
##
BindGlobal("PRIMGRP_InnerGroup",function(inner)
  if inner[1] = "Alt" then
    return AlternatingGroup(inner[2]);
  elif inner[1] = "Sym" then
    return SymmetricGroup(inner[2]);
  elif inner[1] = "PSL" and Length(inner) = 5 then
    return PGOnPointsGroup(inner);
  elif inner[1] = "PSL" then
    return PSL(inner[2],inner[3]);
  elif inner[1] = "PGL" then
    return PGL(inner[2],inner[3]);
  elif inner[1] = "PSigmaL" then
    return PSigmaL(inner[2],inner[3]);
  elif inner[1] = "PGammaL" then
    return PGammaL(inner[2],inner[3]);
  fi;
  Error("PGOnSetsGroup: unknown inner group ",inner[1]);
end);

BindGlobal("PGOnSetsGroup",function(inner,k)
  local g,n,sets,pre,i,c,gens,x,img,pos,j,s,r,prev,h;
  g:=PRIMGRP_InnerGroup(inner);
  n:=LargestMovedPoint(g);
  sets:=Combinations([1..n],k);

  # pre[i][c] counts the k-sets that agree with a given one in their first i-1
  # points and have an i-th point smaller than c
  pre:=List([1..k],i->0*[1..n+1]);
  for i in [1..k] do
    for c in [1..n] do
      pre[i][c+1]:=pre[i][c]+Binomial(n-c,k-i);
    od;
  od;

  gens:=[];
  for x in GeneratorsOfGroup(g) do
    img:=ListPerm(x,n);
    pos:=[];
    for j in [1..Length(sets)] do
      s:=SortedList(img{sets[j]});
      r:=1;
      prev:=0;
      for i in [1..k] do
        r:=r+pre[i][s[i]]-pre[i][prev+1];
        prev:=s[i];
      od;
      pos[j]:=r;
    od;
    Add(gens,PermList(pos));
  od;

  h:=GroupWithGenerators(gens);
  if not IsTransitive(h,[1..Length(sets)]) then
    Error("PGOnSetsGroup: ",inner[1]," is not ",k,"-homogeneous on ",n,
          " points");
  fi;
  return h;
end);

#############################################################################
##
#F  PGPsl( <dim>, <q> ) . . . . the natural projective actions of the L series
#F  PGPgl( <dim>, <q> )
#F  PGPsigmaL( <dim>, <q> )
#F  PGPgammaL( <dim>, <q> )
##
##  PSL, PGL, PSigmaL and PGammaL on the points of PG(dim-1,q), of which there
##  are (q^dim-1)/(q-1).  Unlike PGAlt and PGSym these take arguments, because
##  the degree does not determine the dimension.  They then return functions
##  with argument `deg` and `nr`, which a data file names with its arguments,
##  as ["PSL",dim,q], for PRIMGrp to evaluate on first use.
##
BindGlobal("PGPslOrder",function(dim,q)
  local o,i;
  o:=q^(dim*(dim-1)/2);
  for i in [2..dim] do
    o:=o*(q^i-1);
  od;
  return o/Gcd(dim,q-1);
end);

BindGlobal("PGPsl",function(dim,q)
  local n;
  n:=(q^dim-1)/(q-1);
  return function(deg,nr)
    local t;
    if deg <> n then
      Error("PGPsl(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    t:=2;
    if dim = 2 and q mod 2 = 0 then
      # PSL(2,q) = PGL(2,q) in even characteristic is 3-transitive
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q), 1, "2", [[deg-1,1]], t,
             Concatenation("PSL(",String(dim),",",String(q),")"),
             ["L",[dim,q],1], "psl" ];
  end;
end);

BindGlobal("PGPgl",function(dim,q)
  local n;
  n:=(q^dim-1)/(q-1);
  return function(deg,nr)
    local t;
    if deg <> n then
      Error("PGPgl(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    t:=2;
    if dim = 2 then
      # PGL(2,q) is sharply 3-transitive on the projective line.
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q)*Gcd(dim,q-1), 0, "2", [[deg-1,1]], t,
             Concatenation("PGL(",String(dim),",",String(q),")"),
             ["L",[dim,q],1], "pgl" ];
  end;
end);

BindGlobal("PGPsigmaL",function(dim,q)
  local n,e;
  n:=(q^dim-1)/(q-1);
  e:=Length(Factors(q));
  return function(deg,nr)
    local t;
    if deg <> n then
      Error("PGPsigmaL(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    t:=2;
    if dim = 2 and q mod 2 = 0 then
      # PSL(2,q) = PGL(2,q) in even characteristic is 3-transitive,
      # and PSigmaL contains it
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q)*e, 0, "2", [[deg-1,1]], t,
             Concatenation("PSigmaL(",String(dim),",",String(q),")"),
             ["L",[dim,q],1], "psigmal" ];
  end;
end);

BindGlobal("PGPgammaL",function(dim,q)
  local n,e;
  n:=(q^dim-1)/(q-1);
  e:=Length(Factors(q));
  return function(deg,nr)
    local t;
    if deg <> n then
      Error("PGPgammaL(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    t:=2;
    if dim = 2 then
      # PGL(2,q) is sharply 3-transitive on the projective line,
      # and PGammaL contains it.
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q)*Gcd(dim,q-1)*e, 0, "2", [[deg-1,1]], t,
             Concatenation("PGammaL(",String(dim),",",String(q),")"),
             ["L",[dim,q],1], "pgammal" ];
  end;
end);

#############################################################################
##
#F  PGOnSubspaces( <inner>, <k> ) . . . . . . . . the L series on its k-spaces
#F  PGOnSubspacesGroup( <inner>, <k> )
##
##  The inner group, one of ["PSL",dim,q], ["PGL",dim,q], ["PSigmaL",dim,q] and
##  ["PGammaL",dim,q], acting on the k-dimensional subspaces of GF(q)^dim, for
##  2 <= k <= dim/2.  A k-space and its annihilator give the same permutation
##  group, so a larger k is the action on dim-k, and k = 1 is the inner group's
##  own action on points.  A data file names it ["subspaces", <inner>, k],
##  which is also the entry's field 9.
##
##  The order, the name and the socle are the inner group's, from
##  PRIMGRP_InnerFields.  The suborbits are those of the Grassmann scheme: the
##  k-spaces meeting a given one in a (k-i)-space number
##  q^(i^2) [k,i]_q [dim-k,i]_q, with Gaussian binomials.
##
##  PGOnSubspacesGroup builds the group.  A k-space is the set of the
##  projective points it contains, so the group is built on the points and acts
##  on those sets; the points and the orbit are sorted, as in PGOnSetsGroup.
##
BindGlobal("PRIMGRP_GaussianBinomial",function(n,k,q)
  return Product([0..k-1],i->(q^(n-i)-1)/(q^(i+1)-1));
end);

##  the index of a named projective group over PSL(dim,q)
BindGlobal("PRIMGRP_ProjectiveIndex",function(name,dim,q)
  if name = "PSL" then
    return 1;
  elif name = "PGL" then
    return Gcd(dim,q-1);
  elif name = "PSigmaL" then
    return Length(Factors(q));
  elif name = "PGammaL" then
    return Gcd(dim,q-1)*Length(Factors(q));
  fi;
  Error("unknown projective group ",name);
end);

##  What an inner group gives the entry of an action of it: the number of its
##  points, its order, whether it is simple or solvable, its name and its socle.
BindGlobal("PRIMGRP_InnerFields",function(inner)
  local n,dim,q,idx,flags;
  if inner[1] = "Alt" then
    n:=inner[2];
    return rec(degree:=n, order:=Factorial(n)/2, flags:=1,
               name:=Concatenation("A(",String(n),")"), socle:=["A",n,1]);
  elif inner[1] = "Sym" then
    n:=inner[2];
    return rec(degree:=n, order:=Factorial(n), flags:=0,
               name:=Concatenation("S(",String(n),")"), socle:=["A",n,1]);
  elif inner[1] = "PSL" and Length(inner) = 5 then
    dim:=inner[2];
    q:=inner[3];
    return rec(degree:=(q^dim-1)/(q-1),
               order:=PGPslOrder(dim,q)
                      *Length(PRIMGRP_ProjectiveExtension(dim,q,inner[4])),
               flags:=0, name:=inner[5], socle:=["L",[dim,q],1]);
  elif inner[1] in ["PSL","PGL","PSigmaL","PGammaL"] then
    dim:=inner[2];
    q:=inner[3];
    idx:=PRIMGRP_ProjectiveIndex(inner[1],dim,q);
    if idx = 1 then
      flags:=1;       # simple
    else
      flags:=0;
    fi;
    return rec(degree:=(q^dim-1)/(q-1), order:=PGPslOrder(dim,q)*idx,
               flags:=flags,
               name:=Concatenation(inner[1],"(",String(dim),",",String(q),")"),
               socle:=["L",[dim,q],1]);
  fi;
  Error("unknown inner group ",inner[1]);
end);

BindGlobal("PGOnSubspaces",function(inner,k)
  local dim,q,f;
  dim:=inner[2];
  q:=inner[3];
  f:=PRIMGRP_InnerFields(inner);
  return function(deg,nr)
    Assert(0, 2 <= k and 2*k <= dim
              and deg = PRIMGRP_GaussianBinomial(dim,k,q));
    return [ nr, f.order, f.flags, "2",
             Set(Collected(List([1..k],
               i->q^(i^2)*PRIMGRP_GaussianBinomial(k,i,q)
                         *PRIMGRP_GaussianBinomial(dim-k,i,q)))),
             1, f.name, f.socle, ["subspaces",inner,k] ];
  end;
end);

BindGlobal("PGOnSubspacesGroup",function(inner,k)
  local dim,q,mats,vecs,gens,g,seed,pts;
  dim:=inner[2];
  q:=inner[3];
  Assert(0, 2 <= k and 2*k <= dim);
  if inner[1] in ["PSL","PSigmaL"] then
    mats:=GeneratorsOfGroup(SL(dim,q));
  elif inner[1] in ["PGL","PGammaL"] then
    mats:=GeneratorsOfGroup(GL(dim,q));
  else
    Error("PGOnSubspacesGroup: unknown group ",inner[1]);
  fi;

  vecs:=Set(NormedRowVectors(GF(q)^dim));
  gens:=List(mats,m->Permutation(m,vecs,OnLines));
  if inner[1] in ["PSigmaL","PGammaL"] then
    Add(gens,Permutation(FrobeniusAutomorphism(GF(q)),vecs,
                         function(v,e) return List(v,x->x^e); end));
  fi;
  g:=GroupWithGenerators(gens);

  seed:=Filtered([1..Length(vecs)],i->IsZero(vecs[i]{[k+1..dim]}));
  pts:=Set(Orbit(g,seed,OnSets));
  return Action(g,pts,OnSets);
end);

#############################################################################
##
#F  PGOnSets( <inner>, <k> ) . . . . . . . . . an inner group on its k-subsets
##
##  The inner group, ["Alt",n], ["Sym",n] or a projective one such as
##  ["PSL",dim,q], acting on the k-subsets of its points.  A data file names it
##  ["sets", <inner>, k], which is also the entry's field 9 and is built by
##  PGOnSetsGroup.
##
##  The order, the name and the socle are the inner group's, from
##  PRIMGRP_InnerFields.  The suborbits of Alt(n) and Sym(n) are those of the
##  Johnson scheme.  Those of a projective group are the orbits of the
##  stabiliser of one k-set on the others, computed in the group on the points,
##  which is small.
##
BindGlobal("PRIMGRP_SetsSuborbits",function(inner,k)
  local g,seed,s;
  if inner[1] in ["Alt","Sym"] then
    return PRIMGRP_JohnsonSuborbits(inner[2],k);
  fi;
  g:=PRIMGRP_InnerGroup(inner);
  seed:=[1..k];
  s:=Stabilizer(g,seed,OnSets);
  return Set(Collected(List(Filtered(
           OrbitsDomain(s,Combinations([1..LargestMovedPoint(g)],k),OnSets),
           o->not seed in o),Length)));
end);

BindGlobal("PGOnSets",function(inner,k)
  local f;
  f:=PRIMGRP_InnerFields(inner);
  return function(deg,nr)
    Assert(0, deg = Binomial(f.degree,k));
    return [ nr, f.order, f.flags, "2", PRIMGRP_SetsSuborbits(inner,k), 1,
             f.name, f.socle, ["sets",inner,k] ];
  end;
end);

#############################################################################
##
#F  PGPslExtended( <dim>, <q>, <words>, <name> ) . . between PSL and PGammaL
##
##  The entry for PGOnPointsGroup(["PSL",dim,q,<words>,<name>]).  A data file
##  names it by that list, which is also the entry's field 9.  Every field but
##  the name follows from the words; no rule gives a name like PSL(2,25).2_3,
##  so the description carries it.
##
BindGlobal("PGPslExtended",function(dim,q,words,name)
  local H,t;
  H:=PRIMGRP_ProjectiveExtension(dim,q,words);
  t:=2;
  if dim = 2 and (q mod 2 = 0 or ForAny(H,x->x[1] mod 2 = 1)) then
    # PSL(2,q) for odd q has two orbits on triples, which an element of odd
    # diagonal exponent joins
    t:=3;
  fi;
  return function(deg,nr)
    Assert(0, deg = (q^dim-1)/(q-1));
    return [ nr, PGPslOrder(dim,q)*Length(H), 0, "2", [[deg-1,1]], t, name,
             ["L",[dim,q],1], ["PSL",dim,q,words,name] ];
  end;
end);

#############################################################################
##
#F  PGProductAction4c( <m>, <k>, <els> ) . . . . . . . . .  the product action
##
##  Return a group of O'Nan-Scott type 4c, that is, a subgroup of the full
##  wreath product Sym(m) wreath Sym(k) embedded in Sym(m^k).  A point of the
##  latter is a tuple (c_1,...,c_k) over [1..m], numbered with c_1 most
##  significant as 1 + Sum_i (c_i-1)*m^(k-i).
##
##  The list <els> holds the generators, each as [ p_1, ..., p_k, sigma ],
##  meaning the element that sends (c_1,...,c_k) to (d_1,...,d_k) where
##  d_j = p_i(c_i) for i = j^(sigma^-1).
##
##  Only a minority of the type 4c entries are the whole of a primitive group
##  wreath a transitive one; the rest are proper subgroups, and any of them
##  can be written this way.
##
BindGlobal("PGProductAction4c",function(m,k,els)
  local deg,pw,tuple,index,gens,e,sigma,img,x,c,d,j,i;
  deg:=m^k;
  pw:=List([1..k],i->m^(k-i));
  tuple:=function(x)
    local c,i,r;
    c:=[]; r:=x-1;
    for i in [1..k] do
      c[i]:=QuoInt(r,pw[i])+1;
      r:=r mod pw[i];
    od;
    return c;
  end;
  index:=function(c)
    local i,x;
    x:=1;
    for i in [1..k] do
      x:=x+(c[i]-1)*pw[i];
    od;
    return x;
  end;
  gens:=[];
  for e in els do
    sigma:=e[k+1];
    img:=[];
    for x in [1..deg] do
      c:=tuple(x);
      d:=[];
      for j in [1..k] do
        i:=j^(sigma^-1);
        d[j]:=c[i]^e[i];
      od;
      img[x]:=index(d);
    od;
    Add(gens,PermList(img));
  od;
  return Group(gens);
end);

#############################################################################
##
#F  PRIMGRP_EntryFromDescription( <desc>, <deg>, <nr> ) . . . . an entry, built
##
##  A data file may hold a description of an entry in place of the entry: a
##  list whose first element is a string naming a construction, and whose
##  remaining elements are its arguments.  ["Alt"] is the natural alternating
##  group of the degree it stands at.
##
##  A real entry begins with its own number, so the two are told apart by
##  whether the first element is a string, and nothing has to be marked.
##
##  The names are matched here rather than looked up as globals.  A data file
##  is data: it should be able to ask for one of the constructions the library
##  offers, and for nothing else.
##
BindGlobal("PRIMGRP_EntryFromDescription",function(desc,deg,nr)
  if desc[1] = "Alt" then
    return PGAlt(deg,nr);
  elif desc[1] = "Sym" then
    return PGSym(deg,nr);
  elif desc[1] = "Prime" then
    return PGPrime(desc[2])(deg,nr);
  elif desc[1] = "PSL" and Length(desc) = 5 then
    return PGPslExtended(desc[2],desc[3],desc[4],desc[5])(deg,nr);
  elif desc[1] = "PSL" then
    return PGPsl(desc[2],desc[3])(deg,nr);
  elif desc[1] = "PGL" then
    return PGPgl(desc[2],desc[3])(deg,nr);
  elif desc[1] = "PSigmaL" then
    return PGPsigmaL(desc[2],desc[3])(deg,nr);
  elif desc[1] = "PGammaL" then
    return PGPgammaL(desc[2],desc[3])(deg,nr);
  elif desc[1] = "subspaces" then
    return PGOnSubspaces(desc[2],desc[3])(deg,nr);
  elif desc[1] = "sets" then
    return PGOnSets(desc[2],desc[3])(deg,nr);
  fi;
  Error("unknown construction \"",desc[1],"\" for entry ",nr,
        " of degree ",deg);
end);

#############################################################################
##
##
BindGlobal("PrimGrpLoad",function(deg)
  local s,fname,ind;
  if not IsBound(PRIMGRP[deg]) then
    if not (deg in PRIMRANGE and IsBound(PRIMINDX[deg])) then
      Error("Primitive groups of degree ",deg," are not known!");
    fi;

    ind:=PRIMINDX[deg];
    fname:=Concatenation("gps",String(ind));
    ReadPackage( "primgrp", Concatenation( "data/", fname, ".g" ) );
  fi;
end);

BindGlobal("PRIMGrp",function(deg,nr)
  local l;
  if nr>PRIMLENGTHS[deg] then
    Error("There are only ",PRIMLENGTHS[deg]," groups of degree ",deg,"\n");
  fi;
  PrimGrpLoad(deg);
  l:=PRIMGRP[deg][nr];
  if IsStringRep(l[1]) then
    # The entry is a description of itself rather than itself: a construction
    # named by a string, with its arguments.  Build it, and put the result back
    # so that the next reader finds the entry and not the description.
    l:=PRIMGRP_EntryFromDescription(l,deg,nr);
    PRIMGRP[deg][nr]:=l;
  fi;
  return l;
end);

InstallGlobalFunction(NrPrimitiveGroups, function(deg)
  if deg > 8191 then
    Error("Only groups of degree at most 8191 are known!");
  fi;
  if not IsBound(PRIMLENGTHS[deg]) then
    PrimGrpLoad(deg);
  fi;
  return PRIMLENGTHS[deg];
end);

InstallGlobalFunction(PrimitiveGroupsAvailable,function(deg)
  return deg in PRIMRANGE;
end);

InstallGlobalFunction( PrimitiveGroup, function(deg,num)
local l,g,fac,mats,perms,v,t,filename,strm,r,dim,q,k;

  l:=PRIMGrp(deg,num);

  # special case: Symmetric and Alternating Group
  if l[9]="Alt" then
    g:=AlternatingGroup(deg);
  elif l[9]="Sym" then
    g:=SymmetricGroup(deg);
  elif l[9] in ["psl","pgl","psigmal","pgammal"] then
    # extract the dimension and field from the socle type in field 8
    dim:=l[8][2][1];
    q:=l[8][2][2];
    if l[9] = "psl" then
      g:= PSL(dim,q);
    elif l[9] = "pgl" then
      g:= PGL(dim,q);
    elif l[9] = "psigmal" then
      g:= PSigmaL(dim,q);
    else
      g:= PGammaL(dim,q);
    fi;
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "sets" then
    g:= PGOnSetsGroup(l[9][2], l[9][3]);
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "subspaces" then
    g:= PGOnSubspacesGroup(l[9][2], l[9][3]);
  elif IsList(l[9]) and Length(l[9]) = 5 and l[9][1] = "PSL" then
    g:= PGOnPointsGroup(l[9]);
  elif Length(l[9]) = 2 and l[9][1] = "4c" then
    # product action: the socle width in field 8 gives k, and the degree
    # its k-th root gives m
    k:=l[8][3];
    g:= PGProductAction4c(RootInt(deg,k), k, l[9][2]);
  elif l[4] = "1" then
    # affine type groups described by matrices
    if Length(l[9]) > 0 then
      fac:= Factors(deg);
      mats:=List(l[9],i->ImmutableMatrix(GF(fac[1]),i));
      v:=AsSet(GF(fac[1])^Length(fac));
      perms:=List(mats,i->Permutation(i,v,OnRight));
      t:=First(v,i->not IsZero(i)); # one nonzero translation
                                    #suffices as matrix
                                    # action is irreducible
      Add(perms,Permutation(t,v,function(i,j) return i+j;end));
      g:= Group(perms);
    else
      g:= Image(IsomorphismPermGroup(CyclicGroup(deg)));
    fi;
  else
    # general case: generators given as permutations
    g:= GroupByGenerators( l[9], () );
  fi;

  # now use information from the PRIMGRP entry to prop up the group
  SetPrimitiveIdentification(g,l[1]);
  SetSize(g,l[2]);
  Assert(0, Size(g) = l[2]); # not redundant if g had the size set before
  SetONanScottType(g,l[4]);
  if IsString(l[7]) and Length(l[7])>0 then
    SetName(g,l[7]);
  fi;
  SetSocleTypePrimitiveGroup(g,rec(series:=l[8][1],
                                   parameter:=l[8][2],
                                   width:=l[8][3]));

  if l[3] = 0 then
    SetIsSimpleGroup(g, false);
    SetIsSolvableGroup(g, false);
  elif l[3] = 1 then
    SetIsSimpleGroup(g, true);
    SetIsSolvableGroup(g, false);
  elif l[3] = 2 then
    SetIsSimpleGroup(g, false);
    SetIsSolvableGroup(g, true);
  elif l[3] = 3 then
    SetIsSimpleGroup(g, true);
    SetIsSolvableGroup(g, true);
  fi;

  SetIsAlmostSimpleGroup(g, l[4] = "2");

  SetTransitivity(g, l[6]);
  if deg<=50 then
    SetSimsNo(g, PRIMGRP_SIMSNO[deg-1][num]);
  fi;
  return g;
end );

# local cache for `PrimitiveIdentification':
PRILD:=0;
PGICS:=[];

InstallMethod(PrimitiveIdentification,"generic",true,[IsPermGroup],0,
function(grp)
local dom,deg,PD,s,cand,a,p,s_quot,b,cs,n,beta,alpha,i,ag,bg,q,gl,hom,nr,c,x,conj;
  dom:=MovedPoints(grp);
  if not (IsTransitive(grp,dom) and IsPrimitive(grp,dom)) then
    Error("Group must operate primitively");
  fi;
  deg:=Length(dom);
  # through PRIMGrp, not PRIMGRP[deg]: an entry may be the call that produces
  # it, and PD[i][2] below would be indexing a function.  The upper range has
  # always come this way; the lower one now does too.
  PD:=List([1 .. NrPrimitiveGroups(deg)], t -> PRIMGrp(deg, t));

  if IsNaturalAlternatingGroup(grp) then
    SetSize(grp, Factorial(deg)/2);
  elif IsNaturalSymmetricGroup(grp) then
    SetSize(grp, Factorial(deg));
  fi;

  # size
  s:=Size(grp);
  cand:=Filtered([1..PRIMLENGTHS[deg]],i->PD[i][2]=s);

  #ons
  if Length(cand)>1 and Length(Set(PD{cand},i->i[4]))>1 then
    a:=ONanScottType(grp);
    cand:=Filtered(cand,i->PD[i][4]=a);
  fi;

  # suborbits
  if Length(cand)>1 and Length(Set(PD{cand},i->i[5]))>1 then
    s:=Stabilizer(grp,dom[1]);
    a:=Collected(OrbitLengths(s,dom{[2..Length(dom)]}));
    cand:=Filtered(cand,i->Set(PD[i][5])=Set(a));
  fi;

  # Transitivity
  if Length(cand)>1 and Length(Set(PD{cand},i->i[6]))>1 then
    a:=Transitivity(grp,dom);
    cand:=Filtered(cand,i->PD[i][6]=a);
  fi;

  if Length(cand)>1 then
    # now we need to create the groups
    p:=List(cand,i->PrimitiveGroup(deg,i));

    # in product action case, some tests on the socle quotient.
    if ONanScottType(grp) = "4c" then
     #first we just identify its isomorphism type
      s:= Socle(grp);
      s_quot:= FactorGroup(grp, s);
      a:= IdGroup(s_quot);
      b:= [];
      for i in [1..Length(cand)] do
        b[i]:= IdGroup(FactorGroup(p[i], Socle(p[i])));
      od;
      s:= Filtered([1..Length(cand)], i->b[i] =a);
      cand:= cand{s};
      p:= p{s};
    fi;
  fi;

  # AbelianInvariants
  if Length(cand)>1 then
    a:= AbelianInvariants(grp);
    b:= [];
    for i in [1..Length(cand)] do
      b[i]:= AbelianInvariants(p[i]);
    od;
    s:= Filtered([1..Length(cand)], i->b[i] =a);
    cand:= cand{s};
    p:= p{s};
  fi;

  if Length(cand)>1 then
    # sylow orbits
    gl:=Reversed(PrimeDivisors(Size(grp)));
    while Length(cand)>1 and Length(gl)>0 do
      s:=SylowSubgroup(grp,gl[1]);
      a:=Collected(OrbitLengths(s,MovedPoints(grp)));
      b:=[];
      for i in [1..Length(cand)] do
        s:=SylowSubgroup(p[i],gl[1]);
        b[i]:=Collected(OrbitLengths(s,MovedPoints(p[i])));
      od;
      s:=Filtered([1..Length(cand)],i->b[i]=a);
      cand:=cand{s};
      p:=p{s};
      gl:=gl{[2..Length(gl)]};
    od;
  fi;

  if Length(cand) > 1 then
    # Some further tests for the sylow subgroups
    for q in PrimeDivisors(Size(grp)/Size(Socle(grp))) do
      if q=1 then
        q:=2;
      fi;

      ag:=Image(IsomorphismPcGroup(SylowSubgroup(grp,q)));
      # central series
      a:=List(LowerCentralSeries(ag),Size);
      b:=[];
      for i in [1..Length(cand)] do
        bg:=Image(IsomorphismPcGroup(SylowSubgroup(p[i],q)));
        b[i]:=List(LowerCentralSeries(bg),Size);
      od;
      s:=Filtered([1..Length(cand)],i->b[i]=a);
      cand:=cand{s};
      p:=p{s};

      if Length(cand)>1 then
        # Frattini subgroup
        a:=Size(FrattiniSubgroup(ag));
        b:=[];
        for i in [1..Length(cand)] do
          bg:=Image(IsomorphismPcGroup(SylowSubgroup(p[i],q)));
          b[i]:=Size(FrattiniSubgroup(bg));
        od;
        s:=Filtered([1..Length(cand)],i->b[i]=a);
        cand:=cand{s};
        p:=p{s};
      fi;

      if Length(cand)>1 and Size(ag)<512 then
        # Isomorphism type of 2-Sylow
        a:=IdGroup(ag);
        b:=[];
        for i in [1..Length(cand)] do
          bg:=Image(IsomorphismPcGroup(SylowSubgroup(p[i],q)));
          b[i]:=IdGroup(bg);
        od;
        s:=Filtered([1..Length(cand)],i->b[i]=a);
        cand:=cand{s};
        p:=p{s};
      fi;

    od;
  fi;

  #back for a closer look at the product action groups.
  if Length(cand) > 1 and ONanScottType(grp) = "4c" then
    #just here out of curiosity during testing.
    #Print("cand =", cand, "\n");
    #now we construct the action of the socle quotient as a
    #(necessarily transitive) action on the socle factors.
    s:= Socle(grp);
    cs:= CompositionSeries(s);
    cs:= cs[Length(cs)-1];
    n:= Normalizer(grp, cs);
    beta:= FactorCosetAction(grp, n);
    alpha:= FactorCosetAction(n, ClosureGroup(Centralizer(n, cs), s));
    a:= TransitiveIdentification(Group(KuKGenerators(grp, beta, alpha)));
    b:= [];
    for i in [1..Length(cand)] do
      s:= Socle(p[i]);
      cs:= CompositionSeries(s);
      cs:= cs[Length(cs)-1];
      n:= Normalizer(p[i], cs);
      beta:= FactorCosetAction(p[i], n);
      alpha:= FactorCosetAction(n, ClosureGroup(Centralizer(n, cs), s));
      b[i]:= TransitiveIdentification(Group(KuKGenerators(p[i], beta, alpha)));
    od;
    s:= Filtered([1..Length(cand)], i->b[i]=a);
    cand:= cand{s};
    p:= p{s};
  fi;

  if Length(cand)>1 then
    # Klassen
    a:=Collected(List(ConjugacyClasses(grp:onlysizes),
                      i->[CycleStructurePerm(Representative(i)),Size(i)]));

    # use caching
    if deg<>PRILD then
      PRILD:=deg;
      PGICS:=[];
    fi;

    b:=[];
    for i in [1..Length(cand)] do
      if not IsBound(PGICS[cand[i]]) then
        PGICS[cand[i]]:=Collected(List(ConjugacyClasses(p[i]:onlysizes),
                  j->[CycleStructurePerm(Representative(j)),Size(j)]));
      fi;
      b[i]:=PGICS[cand[i]];
    od;

    s:=Filtered([1..Length(cand)],i->b[i]=a);
    cand:=cand{s};
    p:=p{s};
  fi;

  if Length(cand)>1 and ForAll(p,i->ONanScottType(i)="1")
     and ONanScottType(grp)="1" then
    gl:=Factors(NrMovedPoints(grp));
    gl:=GL(Length(gl),gl[1]);
    hom:=IsomorphismPermGroup(gl);
    s:=List(p,i->Subgroup(gl,LinearActionLayer(i,Pcgs(Socle(i)))));
    b:=Subgroup(gl,LinearActionLayer(grp,Pcgs(Socle(grp))));
    s:=Filtered([1..Length(cand)],
        i->RepresentativeAction(Image(hom,gl),Image(hom,s[i]),Image(hom,b))<>fail);
    cand:=cand{s};
    p:=p{s};
  fi;

  if Length(cand)=1 then
    return cand[1];
  else
    Error("Uh-Oh, this should never happen ",cand);
    return cand[1];
  fi;
end);

InstallMethod(SimsNo,"via `PrimitiveIdentification'",true,[IsPermGroup],0,
function(grp)
local dom;
  dom:=MovedPoints(grp);
  if NrMovedPoints(grp) > 50 then
    Error("SimsNo is defined only for primitive groups of degree <= 50");
  fi;
  if not IsTransitive(grp,dom) and IsPrimitive(grp,dom) then
    Error("Group must operate primitively");
  fi;
  return SimsNo(PrimitiveGroup(Length(dom),PrimitiveIdentification(grp)));
end);

##
#R  IsPrimGrpIterRep
##
DeclareRepresentation("IsPrimGrpIterRep",IsComponentObjectRep,[]);

# function used by the iterator to get the next group or to indicate that
# finished
BindGlobal("PriGroItNext",function(it)
local g;
  it!.next:=fail;
  repeat
    if it!.degi>Length(it!.deg) then
      it!.next:=false;
    else
      g:=PrimitiveGroup(it!.deg[it!.degi],it!.gut[it!.deg[it!.degi]][it!.nr]);
      if ForAll(it!.prop,i->STGSelFunc(i[1](g),i[2])) then
        it!.next:=g;
      fi;
      it!.nr:=it!.nr+1;
      if it!.nr>Length(it!.gut[it!.deg[it!.degi]]) then
        it!.degi:=it!.degi+1;
        it!.nr:=1;
        while it!.degi<=Length(it!.deg) and Length(it!.gut[it!.deg[it!.degi]])=0 do
          it!.degi:=it!.degi+1;
        od;
      fi;
    fi;
  until it!.degi>Length(it!.deg) or it!.next<>fail;
end);

#############################################################################
##
#F  PrimitiveGroupsIterator(arglis,alle)  . . . . . selection function
##
InstallGlobalFunction(PrimitiveGroupsIterator,function(arg)
local arglis,l,deg,pos,mayBeIncomplete,pp,p,requestedDegrees,requestedSizes,
      i,j,a,b,gut,g,grp,nr,RFL,ind,it;
  if Length(arg)=1 and IsList(arg[1]) then
    arglis:=arg[1];
  else
    arglis:=arg;
  fi;
  l:=Length(arglis)/2;
  if not IsInt(l) then
    Error("wrong arguments");
  fi;
  deg:=PRIMRANGE;
  # do we ask for the degree?
  pos:=Filtered([1..l],i->arglis[2*i-1]=NrMovedPoints);
  # The library reaches only up to PRIMRANGE, so a request that may name a
  # degree beyond it cannot be answered in full; that is reported below.
  mayBeIncomplete:= true;
  requestedDegrees:= fail;   # intersection of the degree lists given
  for pp in pos do
    p:=arglis[2*pp];
    if IsInt(p) then
      p:=[p];
    fi;

    if not IsList(p) then
      # a function (wondering, whether anyone will ever use it...)
      deg:= Filtered(deg, p);
      continue;
    fi;

    if requestedDegrees = fail then
      requestedDegrees:= Set(p);
    else
      requestedDegrees:= Intersection(requestedDegrees, p);
    fi;
  od;

  # Only the intersection tells whether the library covers the request:
  # each single list may reach outside PRIMRANGE without a degree being
  # missed, as long as the ones they agree on lie inside.
  if requestedDegrees <> fail then
    mayBeIncomplete:= not IsSubset(PRIMRANGE, requestedDegrees);
    deg:= Intersection(deg, requestedDegrees);
  fi;

  # A primitive group is transitive, so its degree divides its order:
  # order conditions restrict the degree as well, and bound it inside
  # PRIMRANGE as soon as the orders themselves lie there.
  requestedSizes:= fail;
  for ind in [1..l] do
    if arglis[2*ind-1] = Size or arglis[2*ind-1] = Order then
      p:= arglis[2*ind];
      if IsInt(p) then
        p:= [p];
      fi;
      if IsList(p) then
        if requestedSizes = fail then
          requestedSizes:= Set(p);
        else
          requestedSizes:= Intersection(requestedSizes, p);
        fi;
      fi;
    fi;
  od;

  if requestedSizes <> fail then
    mayBeIncomplete:= mayBeIncomplete
                      and not IsSubset(PRIMRANGE, requestedSizes);
    deg:= Filtered(deg, d -> ForAny(requestedSizes, k -> 0 = k mod d));
  elif IsEmpty(pos) then
    Info(InfoWarning,1,"No degree restriction given!\n",
         "#I  A search over the whole library will take a long time!");
  fi;
  gut:=[];
  for i in deg do
    gut[i]:=[1..NrPrimitiveGroups(i)];
  od;

  for i in deg do
    for ind in [1..l] do
      a:=arglis[2*ind-1];
      b:=arglis[2*ind];

      # get all cheap properties first

      if a=NrMovedPoints then
        nr:=0; # done already
      elif a=Size or a=Transitivity or a=ONanScottType then
        if a=Size then
          nr:=2;
        elif a=Transitivity then
          nr:=6;
        elif a=ONanScottType then
          nr:=4;
          if b=1 or b=2 or b=5 then
            b:=String(b);
          elif b=3 then
            b:=["3a","3b"];
          elif b=4 then
            b:=["4a","4b","4c"];
          fi;
        fi;
        gut[i]:=Filtered(gut[i],j->STGSelFunc(PRIMGrp(i,j)[nr],b));
      elif a=IsSimpleGroup or a=IsSimple then
        gut[i]:=Filtered(gut[i],j->STGSelFunc(PRIMGrp(i,j)[3] mod 2=1,b));
      elif a=IsAlmostSimpleGroup or a=IsAlmostSimple then
        # for primitive groups, almost simple means O'Nan-Scott type 2
        gut[i]:=Filtered(gut[i],j->STGSelFunc(PRIMGrp(i,j)[4]="2",b));
      elif a=IsSolvableGroup or a=IsSolvable then
        gut[i]:=Filtered(gut[i],j->STGSelFunc(QuoInt(PRIMGrp(i,j)[3],2)=1,b));
      elif a=SocleTypePrimitiveGroup then
        if IsFunction(b) then
          # for a function we have to translate the list form into records
          RFL:=function(lst)
            return rec(series:=lst[1],parameter:=lst[2],width:=lst[3]);
          end;
          gut[i]:=Filtered(gut[i],j->b(RFL(PRIMGrp(i,j)[8])));
        else
          # otherwise we may bring b into the form we want
          if IsRecord(b) then
            b:=[b];
          fi;
          if IsList(b) and IsRecord(b[1]) then
            b:=List(b,i->[i.series,i.parameter,i.width]);
          fi;
          gut[i]:=Filtered(gut[i],j->PRIMGrp(i,j)[8] in b);
        fi;

      fi;
    od;
  od;

  if mayBeIncomplete then
    Print( "#W  AllPrimitiveGroups: Degree restricted to [ 2 .. ",
           PRIMRANGE[ Length( PRIMRANGE ) ], " ]\n" );
  fi;

  # the rest is hard.

  # find the properties we have not stored
  p:=[];
  for i in [1..l] do
    if not arglis[2*i-1] in
      [NrMovedPoints,Size,Transitivity,ONanScottType,IsSimpleGroup,IsSimple,
       IsAlmostSimpleGroup,IsAlmostSimple,
       IsSolvableGroup,IsSolvable,SocleTypePrimitiveGroup] then
      Add(p,arglis{[2*i-1,2*i]});
    fi;
  od;

  it:=Objectify(NewType(IteratorsFamily,
                        IsIterator and IsPrimGrpIterRep and IsMutable),rec());

  it!.deg:=Immutable(deg);
  i:=1;
  while i<=Length(deg) and Length(gut[deg[i]])=0 do
    i:=i+1;
  od;
  it!.degi:=i;
  it!.nr:=1;
  it!.prop:=MakeImmutable(p);
  it!.gut:=MakeImmutable(gut);
  PriGroItNext(it);
  return it;

end);

InstallMethod(IsDoneIterator,"primitive groups iterator",true,
  [IsPrimGrpIterRep and IsIterator and IsMutable],0,
function(it)
  return it!.next=false or it!.next=fail;
end);

InstallMethod(NextIterator,"primitive groups iterator",true,
  [IsPrimGrpIterRep and IsIterator and IsMutable],0,
function(it)
local g;
  g:=it!.next;
  if g=false or g=fail then
    Error("iterator ran out");
  fi;
  PriGroItNext(it); # next value
  return g;
end);

#############################################################################
##
#F  AllPrimitiveGroups( <fun>, <res>, ... ) . . . . . . . selection function
##
InstallGlobalFunction(AllPrimitiveGroups,function ( arg )
local l,g,it;
  it:=PrimitiveGroupsIterator(arg);
  l:=[];
  for g in it do
    Add(l,g);
  od;
  return l;
end);

#############################################################################
##
#F  OnePrimitiveGroup( <fun>, <res>, ... ) . . . . . . . selection function
##
InstallGlobalFunction(OnePrimitiveGroup,function ( arg )
local l,g,it;
  it:=PrimitiveGroupsIterator(arg);
  if IsDoneIterator(it) then
    return fail;
  else
    return NextIterator(it);
  fi;
end);

# some trivial or useless functions for nitpicking compatibility

BindGlobal("NrAffinePrimitiveGroups",
function(x)
  if x=1 then
    return 1;
  else
   return Length(AllPrimitiveGroups(NrMovedPoints,x,ONanScottType,"1"));
  fi;
end);

BindGlobal("NrSolvableAffinePrimitiveGroups",
  x->Length(AllPrimitiveGroups(NrMovedPoints,x,IsSolvableGroup,true)));

DeclareSynonym("SimsName",Name);

BindGlobal("PrimitiveGroupSims",
function(d,n)
  return OnePrimitiveGroup(NrMovedPoints,d,SimsNo,n);
end);
