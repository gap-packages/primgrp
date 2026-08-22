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
#F  PGAlt( <deg>, <nr> ) . . . . . . . . . . entries we know from the degree
#F  PGSym( <deg>, <nr> )
#F  PGPsl( <deg>, <nr> )
#F  PGPgl( <deg>, <nr> )
##
##  A_n and S_n in their natural action, and PSL(2,q) and PGL(2,q) on the
##  projective line, are determined by the degree, so the data files name them
##  instead of repeating size, suborbits, transitivity, name and socle type
##  for each.  That is ~9300 of the 24558 entries; it also makes it impossible
##  for those fields to drift apart, which is how the socle types of the
##  natural A_n and S_n came to be wrong in degrees >= 2500.
##
##  These build the whole entry, not just the group: PrimitiveGroupsIterator
##  and PrimitiveIdentification read fields 2, 4, 5, 6 and 8 straight out of
##  PRIMGRP without constructing anything.
##
##  Degrees 2, 3 and 4 are excluded -- there A_n and S_n are affine, with an
##  elementary abelian socle -- and their entries are written out in full.
##
BindGlobal("PGAlt",function(deg,nr)
  return [ nr, Factorial(deg)/2, 1, "2", [[deg-1,1]], deg-2,
           Concatenation("Alt(",String(deg),")"), ["A",deg,1], "Alt" ];
end);

BindGlobal("PGSym",function(deg,nr)
  return [ nr, Factorial(deg), 0, "2", [[deg-1,1]], deg,
           Concatenation("Sym(",String(deg),")"), ["A",deg,1], "Sym" ];
end);

#############################################################################
##
#F  PGOnSets( <inner>, <k> ) . . . . . . . a group acting on its own k-subsets
##
##  <inner> names a group in its natural action, and the result is that group
##  acting on the k-subsets of those points.  PrimitiveGroup(2278,3) is
##  Alt(68) on the 2278 = Binomial(68,2) pairs, and PrimitiveGroup(2278,1) is
##  PSL(2,67) on the pairs of the 68 points of the projective line.
##
##  <inner> is one of [ "Alt", n ], [ "Sym", n ], [ "psl", d, q ] or
##  [ "pgl", d, q ].  The inner degree has to be stored: unlike the entries
##  where the action is the natural one, it is not the degree of the entry.
##
##  All of these are k-homogeneous for the k that occur here, so the k-subsets
##  form one orbit and the action is found without enumerating them: the orbit
##  of [1..k] is the whole set.
##
##  The inner group is named by the same kind of tag the data files use for a
##  group in its own right, so anything constructible there can be composed
##  with "take the k-subsets".  PSL(5,2).2 on 2-subsets needs the graph
##  automorphism, which "psl" and "pgl" cannot express but ["cl","L",...] can.
BindGlobal("PRIMGRP_InnerGroup",function(inner)
  if   inner[1] = "Alt" then return AlternatingGroup(inner[2]);
  elif inner[1] = "Sym" then return SymmetricGroup(inner[2]);
  elif inner[1] = "psl" then return PSL(inner[2],inner[3]);
  elif inner[1] = "pgl" then return PGL(inner[2],inner[3]);
  elif inner[1] = "cl"  then
    # These two are defined further down the file; reaching them through
    # ValueGlobal keeps this a forward reference rather than a load-time
    # warning, without moving the definitions apart from their kin.
    return ValueGlobal("PGClassicalGroupOn")(inner[2],inner[3],inner[4],
                                             inner[5],inner[6]);
  elif inner[1] = "clw" then
    return ValueGlobal("PGClassicalWordsOn")(inner[2],inner[3],inner[4],
                                             inner[5],inner[6]);
  fi;
  Error("PGOnSets: unknown inner group ",inner[1]);
end);

#############################################################################
##
#F  PGDiagonal( <tspec>, <ext> ) . . . . . . . . . . . the diagonal action of T
#F  PGDiagonalWords( <tspec>, <words> )
##
##  O'Nan-Scott types 3a and 3b of socle width 2: the socle is T x T acting on
##  the |T| elements of T by g -> x^-1 g y.  The degree is |T|, which is why
##  PrimitiveGroup(7800,6) has socle PSL(2,25) of order 7800.
##
##  The normaliser of that in Sym(|T|) is T^2.(Out(T) x 2), the extra 2 being
##  g -> g^-1, and every entry lies between T^2 and it.  |N/T^2| is 2 to 8
##  here, so an entry is a subgroup of a very small group -- <ext> where the
##  subgroup of that order is unique, <words> in the generators of N where it
##  is not.  Degree 60 has three subgroups of order 2 and needs the latter.
##
BindGlobal("PRIMGRP_DiagonalCache", rec());

BindGlobal("PRIMGRP_SimpleGroupFromSpec",function(tspec)
  if   tspec[1] = "A"    then return AlternatingGroup(tspec[2]);
  elif tspec[1] = "L"    then return PSL(tspec[2],tspec[3]);
  elif tspec[1] = "2A"   then return PSU(tspec[2]+1,tspec[3]);
  elif tspec[1] = "Spor" then
    # only the Mathieu groups occur here, and those are in the GAP library;
    # anything else would need AtlasRep, which primgrp does not depend on
    if tspec[2] = "M(11)" then return MathieuGroup(11); fi;
    if tspec[2] = "M(12)" then return MathieuGroup(12); fi;
    if tspec[2] = "M(22)" then return MathieuGroup(22); fi;
    if tspec[2] = "M(23)" then return MathieuGroup(23); fi;
    if tspec[2] = "M(24)" then return MathieuGroup(24); fi;
  fi;
  Error("no simple group known for ",tspec[1]," ",tspec[2]);
end);

BindGlobal("PGDiagonalNormaliser",function(tspec)
  local key,T,els,act,gens,t,A,a,N,r;
  key:=JoinStringsWithSeparator(List(tspec,String),"_");
  if IsBound(PRIMGRP_DiagonalCache.(key)) then
    return PRIMGRP_DiagonalCache.(key);
  fi;
  T:=PRIMGRP_SimpleGroupFromSpec(tspec);
  els:=AsSSortedList(T);
  # els is sorted, so PositionSorted is a binary search rather than a scan
  act:=f->PermList(List(els,g->PositionSorted(els,f(g))));
  gens:=[];
  for t in GeneratorsOfGroup(T) do
    Add(gens,act(g->g*t));
    Add(gens,act(g->t^-1*g));
  od;
  Add(gens,act(g->g^-1));
  A:=AutomorphismGroup(T);
  for a in GeneratorsOfGroup(A) do
    Add(gens,act(g->Image(a,g)));
  od;
  N:=Group(gens);
  r:=rec(deg:=Length(els), N:=N, T2:=Socle(N));
  PRIMGRP_DiagonalCache.(key):=r;
  return r;
end);

BindGlobal("PGDiagonal",function(tspec,ext)
  local r,hom,Q,cand;
  r:=PGDiagonalNormaliser(tspec);
  if ext=1 then
    return r.T2;
  fi;
  if (Size(r.N)/Size(r.T2)) mod ext <> 0 then
    Error("the diagonal normaliser of ",tspec[1]," realises only ",
          Size(r.N)/Size(r.T2)," outer automorphisms, not ",ext);
  fi;
  hom:=NaturalHomomorphismByNormalSubgroup(r.N,r.T2);
  Q:=ImagesSource(hom);
  cand:=Filtered(List(ConjugacyClassesSubgroups(Q),Representative),
                 s->Size(s)=ext);
  if Length(cand)<>1 then
    Error(Length(cand)," subgroups of order ",ext," above the diagonal socle");
  fi;
  return PreImage(hom,cand[1]);
end);

BindGlobal("PGDiagonalWords",function(tspec,words)
  local r,gens,els;
  r:=PGDiagonalNormaliser(tspec);
  gens:=GeneratorsOfGroup(r.N);
  els:=List(words,w->Product(List(w,e->gens[AbsInt(e)]^SignInt(e)),One(r.N)));
  return ClosureGroup(r.T2,els);
end);

#############################################################################
##
#F  PGAltOnSets( <n>, <k> ) . . . . Alt(n) and Sym(n) on the k-subsets of [1..n]
#F  PGSymOnSets( <n>, <k> )
##
##  Everything about such an entry follows from n and k, so nothing else is
##  stored.  The suborbits are the subdegrees of the Johnson scheme,
##  Binomial(k,i)*Binomial(n-k,i) for i = 1..k, so at n = 45, k = 2 they are
##  2*43 = 86 and 1*903 = 903 -- which is what the file used to spell out
##  beside a 56-digit order.
##
##  The three naming conventions in use, "Alt(45)", "A(5)" and no name at all
##  for 162 of the 362 entries, collapse to the first.
##
BindGlobal("PRIMGRP_JohnsonSuborbits",function(n,k)
  return Set(Collected(List([1..k],i->Binomial(k,i)*Binomial(n-k,i))));
end);

BindGlobal("PGAltOnSets",function(n,k)
  return function(deg,nr)
    Assert(0, deg = Binomial(n,k));
    return [ nr, Factorial(n)/2, 1, "2", PRIMGRP_JohnsonSuborbits(n,k), 1,
             Concatenation("Alt(",String(n),")"), ["A",n,1],
             ["sets",["Alt",n],k] ];
  end;
end);

BindGlobal("PGSymOnSets",function(n,k)
  return function(deg,nr)
    Assert(0, deg = Binomial(n,k));
    return [ nr, Factorial(n), 0, "2", PRIMGRP_JohnsonSuborbits(n,k), 1,
             Concatenation("Sym(",String(n),")"), ["A",n,1],
             ["sets",["Sym",n],k] ];
  end;
end);

BindGlobal("PGOnSetsGroup",function(inner,k)
  local G,pts;
  G:=PRIMGRP_InnerGroup(inner);
  pts:=Orbit(G,[1..k],OnSets);
  if Length(pts) <> Binomial(NrMovedPoints(G),k) then
    Error("PGOnSets: ",inner[1]," is not ",k,"-homogeneous on ",
          NrMovedPoints(G)," points");
  fi;
  return Action(G,pts,OnSets);
end);

#############################################################################
##
#F  PGPrime( <d> ) . . . . . . a subgroup of AGL(1,p) of order p*<d>, p prime
##
##  For prime degree p the affine primitive groups are exactly the subgroups
##  of AGL(1,p) containing the translations, one for each divisor d of p-1.
##  Everything about such an entry follows from p and d, so only d is stored.
##
##  The multiplicative group of GF(p) is cyclic, so its subgroup of order d is
##  unique: any element of order d generates the same group, and the entry is
##  determined by (p,d) rather than by the choice of generator.
##
##  Names follow the convention already used below degree 4096 -- C(p), the
##  dihedral D(2*p), AGL(1,p), and p:d otherwise.  The 7744 entries at degrees
##  4096 to 8191 had no name at all, and gain one here.
##
BindGlobal("PGPrime",function(d)
  return function(deg,nr)
    local a,name,gens,flags,trans;
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
#F  PGPsl( <dim>, <q> ) . . . . . . . . PSL and PGL on projective points
#F  PGPgl( <dim>, <q> )
##
##  PSL(dim,q) and PGL(dim,q) acting on the points of PG(dim-1,q), of which
##  there are (q^dim-1)/(q-1).  Unlike PGAlt and PGSym these take arguments,
##  because the degree does not determine the dimension: 31 points is both the
##  projective line over GF(31)... and PG(4,2).  So the data file says
##  PGPsl(2,13), which yields the function PRIMGrp then applies to (deg,nr).
##
##  Both are 2-transitive on the points; PSL(2,q) is 3-transitive when q is
##  even, where it coincides with PGL(2,q).
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
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q), 1, "2", [[deg-1,1]], t,
             Concatenation("PSL(",String(dim),", ",String(q),")"),
             ["L",[dim,q],1], ["psl",dim,q] ];
  end;
end);

BindGlobal("PGPgl",function(dim,q)
  local n;
  n:=(q^dim-1)/(q-1);
  return function(deg,nr)
    if deg <> n then
      Error("PGPgl(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    return [ nr, PGPslOrder(dim,q)*Gcd(dim,q-1), 0, "2", [[deg-1,1]], 3,
             Concatenation("PGL(",String(dim),", ",String(q),")"),
             ["L",[dim,q],1], ["pgl",dim,q] ];
  end;
end);

#############################################################################
##
#F  PRIMGRP_UnpackMatrices( <p>, <d>, <strings> )
##
##  The matrices of an affine group, packed.  Each d x d matrix over GF(p) is
##  read row by row as a number in base p, that number is written in base 256,
##  and those bytes are base64 encoded.  An 8 x 8 matrix over GF(3) is 101
##  bits, so about 20 characters instead of the 150 its digits take.
##
BindGlobal("PRIMGRP_UnpackMatrices",function(p,d,strings)
  local out, str, bytes, n, c, m, i, j;
  out:=[];
  for str in strings do
    # GAP's Base64String pads the input to a multiple of three bytes and
    # records that in the trailing '=', so the string is kept verbatim.
    bytes:=StringBase64(str);
    n:=0;
    for c in bytes do n:=n*256+INT_CHAR(c); od;
    m:=NullMat(d,d);
    for i in [d,d-1..1] do
      for j in [d,d-1..1] do
        m[i][j]:=n mod p;
        n:=QuoInt(n,p);
      od;
    od;
    Add(out,ImmutableMatrix(GF(p),Z(p)^0*m));
  od;
  return out;
end);

#############################################################################
##
#F  PRIMGRP_AffineVectors( <p>, <d>, <enum> )
##
##  The point set of an affine primitive group of degree p^d, identified with
##  F_p^d.  Two enumerations are in use and they are not the same one:
##
##    "ffe" -- Elements(GF(p)^d), sorted by GAP's order on field elements,
##             which for a prime field runs 0, Z(p)^0, Z(p)^1, ..., that is,
##             by discrete logarithm.  Degrees up to 4095 use this.
##
##    "int" -- by the integer the base-p digits spell out.  Degrees 4096 to
##             8191 use this, since that is how they were imported.
##
##  They agree for p = 2 and p = 3 and differ from p = 5 on, where picking the
##  wrong one silently yields a conjugate of the intended group.  Entries of
##  the second kind therefore store [ "int", <mats> ] rather than bare <mats>.
##
#############################################################################
##
#F  PRIMGRP_AffineAction( <p>, <d>, <enum>, <v>, <mat> )
##
##  The permutation induced on the points by <mat>.
##
##  For the "int" enumeration the base-p digits of a point number are the
##  vector itself, so the image number is read off the image vector and no
##  search is needed.  Permutation(mat,v,OnRight) instead looks each image up
##  in <v>, which is a linear scan over a hand-built list: 4985 ms against
##  17 ms for one generator of PrimitiveGroup(6561,1).
##
BindGlobal("PRIMGRP_AffineTranslation",function(p,d,enum,v,t)
  local pw,img,i,w,x,j;
  if enum <> "int" then
    return Permutation(t,v,function(a,b) return a+b; end);
  fi;
  pw:=List([1..d],j->p^(d-j));
  img:=EmptyPlist(p^d);
  for i in [1..p^d] do
    w:=v[i]+t;
    x:=1;
    for j in [1..d] do x:=x+IntFFE(w[j])*pw[j]; od;
    img[i]:=x;
  od;
  return PermList(img);
end);

BindGlobal("PRIMGRP_AffineAction",function(p,d,enum,v,mat)
  local pw,img,i,w,x,j;
  if enum <> "int" then
    return Permutation(mat,v,OnRight);
  fi;
  pw:=List([1..d],j->p^(d-j));
  img:=EmptyPlist(p^d);
  for i in [1..p^d] do
    w:=v[i]*mat;
    x:=1;
    for j in [1..d] do x:=x+IntFFE(w[j])*pw[j]; od;
    img[i]:=x;
  od;
  return PermList(img);
end);

BindGlobal("PRIMGRP_AffineVectors",function(p,d,enum)
  local out,i,x,w,j;
  if enum = "ffe" then
    return Elements(GF(p)^d);
  elif enum <> "int" then
    Error("unknown affine point enumeration \"",enum,"\"");
  fi;
  out:=[];
  for i in [0..p^d-1] do
    x:=i; w:=[];
    for j in [1..d] do
      Add(w,x mod p);
      x:=QuoInt(x,p);
    od;
    Add(out,Z(p)^0*Reversed(w));
  od;
  return out;
end);

#############################################################################
##
#F  PGPsigmaL( <dim>, <q> ) . . . PSigmaL and PGammaL on projective points
#F  PGPgammaL( <dim>, <q> )
##
##  The extensions of PSL(dim,q) by field automorphisms, and by field and
##  diagonal automorphisms, on the same (q^dim-1)/(q-1) points.  For q = p^e,
##
##      |PSigmaL(dim,q)| = |PSL(dim,q)| * e
##      |PGammaL(dim,q)| = |PGL(dim,q)| * e
##
##  and both are 2-transitive on the points, PGammaL 3-transitive when dim is
##  2 because it then contains the sharply 3-transitive PGL(2,q).  Read off
##  the 25 PSigmaL and 43 PGammaL entries the library already names, which the
##  formulas reproduce exactly.
##
##  These apply only when the degree really is (q^dim-1)/(q-1).  A name is not
##  enough to go on: PrimitiveGroup(496,7) is called "PGammaL(2, 32)" and acts
##  on 496 points, not on the 33 of the projective line, and is not even
##  2-transitive there.
##
#############################################################################
##
#F  PGProductAction( <m>, <b>, <topgens> )
##
##  A group of O'Nan-Scott type 4c: PrimitiveGroup(m,b) wreath <top> in the
##  product action on m^k points, where <top> is the transitive group of
##  degree k generated by <topgens>.
##
##  The top group is written out rather than named, because naming it would
##  mean calling TransitiveGroup and so make transgrp a dependency of the data
##  itself.  It has degree at most 5 here, so its generators cost nothing.
##
##  This replaces only the generators.  The rest of the entry stays as it is,
##  so PrimitiveGroupsIterator still answers from the file without building
##  anything -- and the generators are all but the whole of such an entry.
##
#############################################################################
##
#F  PGProductAction4c( <m>, <k>, <els> )
##
##  A group of O'Nan-Scott type 4c given by its own generators, written as
##  elements of Sym(m) wreath Sym(k) rather than as permutations of the m^k
##  points.  Each element of <els> is [ p_1, ..., p_k, sigma ]: it sends the
##  point with coordinates (c_1,...,c_k) to the one with coordinates
##  d_j = p_i(c_i) where i = j^(sigma^-1).
##
##  Points are numbered with c_1 most significant, so point
##  1 + Sum_i (c_i-1)*m^(k-i) is (c_1,...,c_k).  Only 160 type 4c entries are
##  the whole of PrimitiveGroup(m,b) wreath a transitive group; the rest are
##  proper subgroups, and this writes down any of them.
##
BindGlobal("PGProductAction4c",function(m,k,els)
  local deg,pw,tuple,index,gens,e,sigma,img,x,c,d,j,i;
  deg:=m^k;
  pw:=List([1..k],i->m^(k-i));
  tuple:=function(x)
    local c,i,r;
    c:=[]; r:=x-1;
    for i in [1..k] do c[i]:=QuoInt(r,pw[i])+1; r:=r mod pw[i]; od;
    return c;
  end;
  index:=function(c)
    local i,x;
    x:=1;
    for i in [1..k] do x:=x+(c[i]-1)*pw[i]; od;
    return x;
  end;
  gens:=[];
  for e in els do
    sigma:=e[k+1];
    img:=[];
    for x in [1..deg] do
      c:=tuple(x); d:=[];
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

BindGlobal("PGProductActionGroup",function(m,b,topgens)
  return WreathProductProductAction(PrimitiveGroup(m,b), Group(topgens));
end);

#############################################################################
##
#F  PGClassicalGroup( <ser>, <d>, <q>, <ext> )
##
##  An almost simple group of O'Nan-Scott type 2 whose socle is the classical
##  group <ser>_<d>(<q>), acting on the natural geometry: the projective
##  points for series L and C, the isotropic points for 2A, the singular
##  points for B, D and 2D.
##
##  Only the big group is built.  GAP's Omega(0,7,3) and GO(0,7,3) do not use
##  the same Gram matrix, so their singular points are different sets and no
##  one action carries both; taking the socle of the big group instead makes
##  that mismatch impossible.
##
##  Everything between the socle T and A = N(T) is primitive of the same
##  degree, so the subgroups of A/T of order <ext> are exactly the library
##  entries of that order.  <ext> is recorded only where that subgroup is
##  unique, which is what makes the entry recoverable.
##
BindGlobal("PRIMGRP_ClassicalCache", rec());

##  <variety> picks which of the group's varieties carries the action:
##
##    "points" -- the projective, isotropic or singular 1-spaces.
##    "lines"  -- the totally singular 2-spaces.
##
##  For B(2,q) both have (q^2+1)(q+1) elements and the two actions are
##  inequivalent for odd q, the quadrangle Q(4,q) being self-dual only in even
##  characteristic.  That is why degrees 40, 156, 400, 820, 1464, 2380, 5220
##  and 7240 each carry two entries of equal order, equal socle and equal
##  suborbits.  They differ in the permutation character, and in the
##  isomorphism type of a point stabiliser -- the two maximal parabolics.
##
BindGlobal("PGClassicalNormaliserOn",function(ser,d,q,variety)
  local key,n,Q,M,kind,all,pts,m,bil,u,cand,L,p,frob,act,onset,A,r,gens;
  key:=Concatenation(ser,"_",String(d),"_",String(q),"_",variety);
  if IsBound(PRIMGRP_ClassicalCache.(key)) then
    return PRIMGRP_ClassicalCache.(key);
  fi;
  if   ser="L"  then n:=d;     Q:=q;   M:=GL(d,q);         kind:="all";
  elif ser="C"  then n:=2*d;   Q:=q;   M:=Sp(2*d,q);       kind:="all";
  elif ser="2A" then n:=d+1;   Q:=q^2; M:=GU(d+1,q);       kind:="sesq";
  elif ser="B"  then n:=2*d+1; Q:=q;   M:=GO(0,2*d+1,q);   kind:="quad";
  elif ser="D"  then n:=2*d;   Q:=q;   M:=GO(1,2*d,q);     kind:="quad";
  elif ser="2D" then n:=2*d;   Q:=q;   M:=GO(-1,2*d,q);    kind:="quad";
  else Error("no natural action known for series ",ser);
  fi;
  all:=NormedRowVectors(GF(Q)^n);
  if kind="all" then
    pts:=all;
  elif kind="sesq" then
    m:=InvariantSesquilinearForm(M).matrix;
    pts:=Filtered(all,v->v*m*List(v,x->x^RootInt(Q))=0*Z(Q)^0);
  else
    m:=InvariantQuadraticForm(M).matrix;
    pts:=Filtered(all,v->v*m*v=0*Z(Q)^0);
  fi;
  p:=SmallestRootInt(Q);
  if variety="points" then
    onset:=OnLines;
    act:=function(v,e) return OnLines(List(v,x->x^e),1); end;
  elif variety="lines" then
    if kind="all" then
      Error("no singular lines for series ",ser);
    fi;
    # transitive on totally singular lines, so one orbit is the whole variety
    bil:=m+TransposedMat(m);
    u:=pts[1];
    cand:=Filtered(pts,v->v*bil*u=0*Z(Q)^0 and RankMat([u,v])=2);
    if IsEmpty(cand) then
      Error("no totally singular line through the first point of ",ser);
    fi;
    L:=TriangulizedMat([u,cand[1]]);
    pts:=Orbit(M,L,OnSubspacesByCanonicalBasis);
    onset:=OnSubspacesByCanonicalBasis;
    act:=function(L,e) return TriangulizedMat(List(L,x->List(x,y->y^e))); end;
  else
    Error("unknown variety \"",variety,"\"");
  fi;
  gens:=GeneratorsOfGroup(Action(M,pts,onset));
  # Frobenius permutes the points only when the Gram matrix is fixed by it,
  # which fails for GO(-1,8,4) and its kin over non-prime fields of even
  # characteristic.  Leave it out there: A comes out smaller, and asking for
  # an <ext> it can no longer realise is refused below rather than guessed.
  frob:=Permutation(p,pts,act);
  if frob<>fail then
    gens:=Concatenation(gens,[frob]);
  fi;
  A:=Group(gens);
  r:=rec(deg:=Length(pts), A:=A, T:=Socle(A));
  PRIMGRP_ClassicalCache.(key):=r;
  return r;
end);

#############################################################################
##
#F  PGClassicalWords( <ser>, <d>, <q>, <words> )
##
##  The same groups as PGClassicalGroup, for the degrees where the order does
##  not single one out: PSigmaL(2,25) and PSL(2,25).2_3 are both index 2 over
##  PSL(2,25) on 26 points, and share socle, order and suborbits.
##
##  <words> name the generators to adjoin to the socle, as words in the
##  generators of A: a list of integers, negative for an inverse.  Naming them
##  rather than searching keeps reading cheap -- picking the right one out of
##  the candidates needs an isomorphism invariant, which is fine once during
##  conversion but not on every call.
##
BindGlobal("PGClassicalWordsOn",function(ser,d,q,words,variety)
  local r,gens,els;
  r:=PGClassicalNormaliserOn(ser,d,q,variety);
  gens:=GeneratorsOfGroup(r.A);
  els:=List(words,w->Product(List(w,e->gens[AbsInt(e)]^SignInt(e)),One(r.A)));
  return ClosureGroup(r.T,els);
end);

BindGlobal("PGClassicalWords",function(ser,d,q,words)
  return PGClassicalWordsOn(ser,d,q,words,"points");
end);

BindGlobal("PGClassicalNormaliser",function(ser,d,q)
  return PGClassicalNormaliserOn(ser,d,q,"points");
end);

BindGlobal("PGClassicalGroupOn",function(ser,d,q,ext,variety)
  local r,hom,Q,cand;
  r:=PGClassicalNormaliserOn(ser,d,q,variety);
  if ext=1 then
    return r.T;
  fi;
  if Size(r.A)/Size(r.T) mod ext <> 0 then
    Error(ser,"(",d,",",q,") realises only ",Size(r.A)/Size(r.T),
          " outer automorphisms on its points, not ",ext);
  fi;
  hom:=NaturalHomomorphismByNormalSubgroup(r.A,r.T);
  Q:=ImagesSource(hom);
  cand:=Filtered(List(ConjugacyClassesSubgroups(Q),Representative),
                 s->Size(s)=ext);
  if Length(cand)<>1 then
    Error(Length(cand)," subgroups of order ",ext," in the outer group of ",
          ser,"(",d,",",q,")");
  fi;
  return PreImage(hom,cand[1]);
end);

BindGlobal("PGClassicalGroup",function(ser,d,q,ext)
  return PGClassicalGroupOn(ser,d,q,ext,"points");
end);

BindGlobal("PGPsigmaL",function(dim,q)
  local n,e;
  n:=(q^dim-1)/(q-1);
  e:=Length(Factors(q));
  return function(deg,nr)
    if deg <> n then
      Error("PGPsigmaL(",dim,",",q,") describes degree ",n,", not ",deg);
    fi;
    return [ nr, PGPslOrder(dim,q)*e, 0, "2", [[deg-1,1]], 2,
             Concatenation("PSigmaL(",String(dim),", ",String(q),")"),
             ["L",[dim,q],1], ["psigmal",dim,q] ];
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
      t:=3;
    fi;
    return [ nr, PGPslOrder(dim,q)*Gcd(dim,q-1)*e, 0, "2", [[deg-1,1]], t,
             Concatenation("PGammaL(",String(dim),", ",String(q),")"),
             ["L",[dim,q],1], ["pgammal",dim,q] ];
  end;
end);

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
  local r;
  if nr>PRIMLENGTHS[deg] then
    Error("There are only ",PRIMLENGTHS[deg]," groups of degree ",deg,"\n");
  fi;
  PrimGrpLoad(deg);
  # An entry may be a function rather than a list, for the groups the degree
  # alone determines: the data file says PGSym, not a spelled-out entry whose
  # size field is Factorial(8191).  Evaluate on demand, and keep the result.
  r:=PRIMGRP[deg][nr];
  if IsFunction(r) then
    r:=r(deg,nr);
    PRIMGRP[deg][nr]:=r;
  fi;
  return r;
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
  return deg <= 8191;
end);

InstallGlobalFunction( PrimitiveGroup, function(deg,num)
local l,g,gens,enum,fac,mats,perms,v,t;

  l:=PRIMGrp(deg,num);

  # special case: Symmetric and Alternating Group
  if l[9]="Alt" then
    g:=AlternatingGroup(deg);
    SetName(g,Concatenation("A(",String(deg),")"));
  elif l[9]="Sym" then
    g:=SymmetricGroup(deg);
    SetName(g,Concatenation("S(",String(deg),")"));
  elif l[9] = "psl" then
    g:= PSL(2, deg-1);
    SetName(g, Concatenation("PSL(2,", String(deg-1),")"));
  elif l[9] = "pgl" then
    g:= PGL(2, deg-1);
    SetName(g, Concatenation("PGL(2,", String(deg-1), ")"));
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "psl" then
    g:= PSL(l[9][2], l[9][3]);
    SetName(g, Concatenation("PSL(", String(l[9][2]), ",", String(l[9][3]), ")"));
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "pgl" then
    g:= PGL(l[9][2], l[9][3]);
    SetName(g, Concatenation("PGL(", String(l[9][2]), ",", String(l[9][3]), ")"));
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "psigmal" then
    g:= PSigmaL(l[9][2], l[9][3]);
    SetName(g, Concatenation("PSigmaL(", String(l[9][2]), ",", String(l[9][3]), ")"));
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "diag" then
    g:= PGDiagonal(l[9][2], l[9][3]);
    if IsString(l[7]) and Length(l[7])>0 then SetName(g,l[7]); fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "diagw" then
    g:= PGDiagonalWords(l[9][2], l[9][3]);
    if IsString(l[7]) and Length(l[7])>0 then SetName(g,l[7]); fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "sets" then
    g:= PGOnSetsGroup(l[9][2], l[9][3]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 6 and l[9][1] = "clw" then
    g:= PGClassicalWordsOn(l[9][2], l[9][3], l[9][4], l[9][5], l[9][6]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 6 and l[9][1] = "cl" then
    g:= PGClassicalGroupOn(l[9][2], l[9][3], l[9][4], l[9][5], l[9][6]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 5 and l[9][1] = "clw" then
    g:= PGClassicalWords(l[9][2], l[9][3], l[9][4], l[9][5]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 5 and l[9][1] = "cl" then
    g:= PGClassicalGroup(l[9][2], l[9][3], l[9][4], l[9][5]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 4 and l[9][1] = "4c" then
    g:= PGProductAction4c(l[9][2], l[9][3], l[9][4]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 4 and l[9][1] = "pa" then
    g:= PGProductActionGroup(l[9][2], l[9][3], l[9][4]);
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    fi;
    SetSize(g,l[2]);
  elif IsList(l[9]) and Length(l[9]) = 3 and l[9][1] = "pgammal" then
    g:= PGammaL(l[9][2], l[9][3]);
    SetName(g, Concatenation("PGammaL(", String(l[9][2]), ",", String(l[9][3]), ")"));
  elif l[4] = "1" then
    gens:=l[9];
    enum:="ffe";
    fac:=Factors(deg);
    if Length(gens) = 3 and IsString(gens[1]) and gens[1] = "b64" then
      enum:=gens[2];
      gens:=PRIMGRP_UnpackMatrices(fac[1],Length(fac),gens[3]);
    elif Length(gens) = 2 and IsString(gens[1]) then
      enum:=gens[1];             # see PRIMGRP_AffineVectors below
      gens:=gens[2];
    fi;
    if Length(gens) > 0 then
      fac:= Factors(deg);
      mats:=List(gens,i->ImmutableMatrix(GF(fac[1]),i));
      v:=PRIMGRP_AffineVectors(fac[1],Length(fac),enum);
      perms:=List(mats,i->PRIMGRP_AffineAction(fac[1],Length(fac),enum,v,i));
      t:=First(v,i->not IsZero(i)); # one nonzero translation
                                    #suffices as matrix
                                    # action is irreducible
      Add(perms,PRIMGRP_AffineTranslation(fac[1],Length(fac),enum,v,t));
      g:= Group(perms);
      SetSize(g, l[2]);
    else
      # No matrix generators: the group is the translations alone, which is
      # regular, and a regular elementary abelian group is primitive only in
      # dimension one -- so deg is prime and this is the cyclic group of
      # order deg in its regular action.  All 1026 entries that reach here
      # are of that shape.  Say so directly: IsomorphismPermGroup(
      # CyclicGroup(deg)) took 500 ms at degree 8191 and guaranteed nothing
      # about the representation it returned, and CyclicGroup(IsPermGroup,
      # deg) leaves the degree to chance.
      g:= Group(CycleFromList([1..deg]));
      SetSize(g, deg);
    fi;
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g, l[7]);
    fi;
  else
    g:= GroupByGenerators( l[9], () );
    if IsString(l[7]) and Length(l[7])>0 then
      SetName(g,l[7]);
    #else
    #  SetName(g,Concatenation("p",String(deg),"n",String(num)));
    fi;
    SetSize(g,l[2]);
  fi;
  SetPrimitiveIdentification(g,l[1]);
  SetONanScottType(g,l[4]);
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
  SetTransitivity(g, l[6]);
  if deg<=50 then
    SetSimsNo(g,PRIMGRP_SIMSNO[deg-1][num]);
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
  # via PRIMGrp, not PRIMGRP directly: entries may still be unevaluated
  PD:=List([1 .. NrPrimitiveGroups(deg)], t -> PRIMGrp(deg, t));

  if IsNaturalAlternatingGroup(grp) then
    SetSize(grp, Factorial(deg)/2);
  elif IsNaturalSymmetricGroup(grp) then
    SetSize(grp, Factorial(deg));
  fi;

  s:=Size(grp);

  # size
  cand:=Filtered([1..PRIMLENGTHS[deg]],i->PD[i][2]=s);

  #ons
  if Length(cand)>1 and Length(Set(PD{cand},i->i[4]))>1 then
    a:=ONanScottType(grp);
    cand:=Filtered(cand,i->PD[i][4]=a);
  fi;

  # suborbits
  if Length(cand)>1 and Length(Set(PD{cand},i->i[5]))>1 then
    a:=Collected(List(Orbits(Stabilizer(grp,dom[1]),dom{[2..Length(dom)]}),
                      Length));
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
    gl:=Reversed(Set(Factors(Size(grp))));
    while Length(cand)>1 and Length(gl)>0 do
      a:=Collected(List(Orbits(SylowSubgroup(grp,gl[1]),MovedPoints(grp)),
                        Length));
      b:=[];
      for i in [1..Length(cand)] do
        b[i]:=Collected(List(Orbits(SylowSubgroup(p[i],gl[1]),
                                    MovedPoints(p[i])),
                          Length));
      od;
      s:=Filtered([1..Length(cand)],i->b[i]=a);
      cand:=cand{s};
      p:=p{s};
      gl:=gl{[2..Length(gl)]};
    od;
  fi;

  if Length(cand) > 1 then
    # Some further tests for the sylow subgroups
    for q in Set(Factors(Size(grp)/Size(Socle(grp)))) do
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
local arglis,i,j,a,b,l,p,deg,gut,g,grp,nr,f,RFL,ind,it;
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
  p:=Position(arglis,NrMovedPoints);
  if p<>fail then
    p:=arglis[p+1];
    if IsInt(p) then
      f:=not p in deg;
      p:=[p];
    fi;
    if IsList(p) then
      f:=not IsSubset(deg,Difference(p,[1]));
      deg:=Intersection(deg,p);
    else
      # b is a function (wondering, whether anyone will ever use it...)
      f:=true;
      deg:=Filtered(deg,p);
    fi;
  else
    f:=true; #warnung weil kein Degree angegeben ?
    b:=true;
    for a in [Size,Order] do
      p:=Position(arglis,a);
      if p<>fail then
        p:=arglis[p+1];
        if IsInt(p) then
          p:=[p];
        fi;

        if IsList(p) then
          deg := Filtered( deg,
               d -> ForAny( p, k -> 0 = k mod d ) );
          b := false;
          f := not IsSubset( PRIMRANGE, p );
        fi;
      fi;
    od;
    if b then
      Info(InfoWarning,1,"No degree restriction given!\n",
           "#I  A search over the whole library will take a long time!");
    fi;
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

  if f then
    Print( "#W  AllPrimitiveGroups: Degree restricted to [ 2 .. ",
           PRIMRANGE[ Length( PRIMRANGE ) ], " ]\n" );
  fi;

  # the rest is hard.

  # find the properties we have not stored
  p:=[];
  for i in [1..l] do
    if not arglis[2*i-1] in
      [NrMovedPoints,Size,Transitivity,ONanScottType,IsSimpleGroup,IsSimple,
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
