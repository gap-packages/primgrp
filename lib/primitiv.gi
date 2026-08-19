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
      perms:=List(mats,i->Permutation(i,v,OnRight));
      t:=First(v,i->not IsZero(i)); # one nonzero translation
                                    #suffices as matrix
                                    # action is irreducible
      Add(perms,Permutation(t,v,function(i,j) return i+j;end));
      g:= Group(perms);
      SetSize(g, l[2]);
    else
      g:= Image(IsomorphismPermGroup(CyclicGroup(deg)));
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
