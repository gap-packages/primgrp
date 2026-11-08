#############################################################################
##
#W  constructions.gi  GAP primitive groups library          Alexander Hulpke
##
##
#Y  Copyright (C)  1999, School Math.&Comp. Sci., University of St Andrews
##
##  The constructions a data file may name in place of an entry.
##

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
BindGlobal("PGPrime",function(d)
  return function(deg,nr)
    local name,gens,flags,trans;
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
      flags := 3;  # simple and solvable
    else
      flags := 2;  # solvable
      gens:=[ [ [ Z(deg)^((deg-1)/d) ] ] ];
      if d = 2 then
        name:=Concatenation("D(2*",String(deg),")");
      elif d = deg-1 then
        name:=Concatenation("AGL(1, ",String(deg),")");
      else
        name:=Concatenation(String(deg),":",String(d));
      fi;
    fi;
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
#F  PRIMGRP_ProjectiveExtension( <dim>, <q>, <auts> ) . . . their closure
##
##  Write q = p^f, and g = gcd(dim,q-1).  In PGammaL(dim,q) the diagonal
##  matrix delta = diag(Z(q),1,...,1) and the Frobenius phi give coset
##  representatives of PSL(dim,q): every coset is the one of delta^i phi^j for
##  integers 0 <= i < g and 0 <= j < f, and the pair [i,j] names that coset.
##  Nothing outside those ranges is needed, delta^g and phi^f lying in
##  PSL(dim,q).
##
##  Composing representatives, first delta^i1 phi^j1 and then delta^i2 phi^j2,
##  lands in the coset named by [i1 + i2*p^(f-j1) mod g, j1 + j2 mod f].  A
##  group G with PSL(dim,q) <= G <= PGammaL(dim,q) is a union of cosets closed
##  under that, and <auts> is a list of pairs naming some of them.
##
##  Return the list of pairs naming the cosets of the smallest such G, that is
##  of the group the pairs in <auts> generate.  Their number is the index of
##  PSL(dim,q) in G, so G has that many times the order of PSL(dim,q).
##
BindGlobal("PRIMGRP_ProjectiveExtension",function(dim,q,auts)
  local p,f,g,mul,H,queue,x,y,z;
  p:=Characteristic(GF(q));
  f:=Length(Factors(q));
  g:=Gcd(dim,q-1);
  mul:=function(a,b)
    return [ (a[1]+b[1]*PowerModInt(p,(f-a[2]) mod f,g)) mod g,
             (a[2]+b[2]) mod f ];
  end;
  H:=[[0,0]];
  queue:=[[0,0]];
  while queue <> [] do
    x:=Remove(queue);
    for y in auts do
      z:=mul(x,[y[1] mod g,y[2] mod f]);
      if not z in H then
        Add(H,z);
        Add(queue,z);
      fi;
    od;
  od;
  return Set(H);
end);

#############################################################################
##
#F  PGPslExtendedGroup( <inner> ) . PSL extended by automorphisms, on points
##
##  <inner> is ["PSL",dim,q,<auts>,<name>]: PSL(dim,q) on the points of
##  PG(dim-1,q), extended by one element for each automorphism in <auts>,
##  diag(Z(q)^i,1,...,1) followed by the j-th power of the Frobenius for the
##  pair [i,j].  Every group with PSL(dim,q) <= G <= PGammaL(dim,q) arises
##  this way, those that are none of the four named ones included.  The points
##  are the sorted normed vectors.
##
BindGlobal("PGPslExtendedGroup",function(inner)
  local dim,q,p,vecs,gens,a,e,u,one;
  dim:=inner[2];
  q:=inner[3];
  p:=Characteristic(GF(q));
  vecs:=Set(NormedRowVectors(GF(q)^dim));
  gens:=List(GeneratorsOfGroup(SL(dim,q)),m->Permutation(m,vecs,OnLines));
  one:=One(GF(q));
  for a in inner[4] do
    # delta^i phi^j raises every coordinate to the power e = p^j and multiplies
    # the first by Z(q)^(i*e).  A normed vector stays normed while its first
    # coordinate is zero; otherwise that coordinate becomes Z(q)^(i*e), and
    # dividing the vector by it norms the vector again.
    e:=p^a[2];
    u:=Z(q)^(a[1]*e);
    Add(gens,Permutation((),vecs,function(v,x)
      local w;
      w:=List(v,c->c^e);
      if IsZero(w[1]) then
        return w;
      fi;
      w:=w/u;
      w[1]:=one;
      return w;
    end));
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
##  The points are all k-subsets of [1..n], in the order Combinations gives
##  them.  That order is fixed; one given by Orbit would not be, and which
##  permutation group comes out depends on it.  The result is transitive only
##  if the group is k-homogeneous, so that is checked.
##
BindGlobal("PRIMGRP_InnerGroup",function(inner)
  if inner[1] = "Alt" then
    return AlternatingGroup(inner[2]);
  elif inner[1] = "Sym" then
    return SymmetricGroup(inner[2]);
  elif inner[1] = "PSL" and Length(inner) = 5 then
    return PGPslExtendedGroup(inner);
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
  local g,sets,h;
  g:=PRIMGRP_InnerGroup(inner);
  sets:=Combinations([1..LargestMovedPoint(g)],k);
  h:=Action(g,sets,OnSets);
  if not IsTransitive(h,[1..Length(sets)]) then
    Error("PGOnSetsGroup: ",inner[1]," is not ",k,"-homogeneous on ",
          LargestMovedPoint(g)," points");
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
#F  PGPslExtended( <dim>, <q>, <auts>, <name> ) . .  between PSL and PGammaL
##
##  The entry for PGPslExtendedGroup(["PSL",dim,q,<auts>,<name>]).  A data
##  file names it by that list, which is also the entry's field 9.  Every
##  field but the name follows from <auts>; no rule gives a name like
##  PSL(2,25).2_3, so the description carries it.
##
BindGlobal("PGPslExtended",function(dim,q,auts,name)
  local H,t;
  H:=PRIMGRP_ProjectiveExtension(dim,q,auts);
  t:=2;
  if dim = 2 and (q mod 2 = 0 or ForAny(H,x->x[1] mod 2 = 1)) then
    # PSL(2,q) for odd q has two orbits on triples, which an element of odd
    # diagonal exponent joins
    t:=3;
  fi;
  return function(deg,nr)
    Assert(0, deg = (q^dim-1)/(q-1));
    return [ nr, PGPslOrder(dim,q)*Length(H), 0, "2", [[deg-1,1]], t, name,
             ["L",[dim,q],1], ["PSL",dim,q,auts,name] ];
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
