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
##  An entry that is a function is evaluated by PRIMGrp on first use. This
##  keeps Factorial(deg) from being computed for every degree in a data file
##  merely because the file was read.
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
#F  PGPsl( <dim>, <q> ) . . . . the natural projective actions of the L series
#F  PGPgl( <dim>, <q> )
#F  PGPsigmaL( <dim>, <q> )
#F  PGPgammaL( <dim>, <q> )
##
##  PSL, PGL, PSigmaL and PGammaL on the points of PG(dim-1,q), of which there
##  are (q^dim-1)/(q-1).  Unlike PGAlt and PGSym these take arguments, because
##  the degree does not determine the dimension. They then return functions
##  with argument `deg` and `nr` that are inserted as entries into PRIMGRP.
##  These functions are then invoked as needed by PRIMGrp.
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
##
BindGlobal("PrimGrpLoad",function(deg)
  local s,fname,ind;
  if not IsBound(PRIMGRP[deg]) then
    if deg > 4095 then
      Error("This method is not for primitive groups of degree greater than 4095!");
    fi;
    if not (deg in PRIMRANGE and IsBound(PRIMINDX[deg])) then
      Error("Primitive groups of degree ",deg," are not known!");
    fi;

    ind:=PRIMINDX[deg];
    fname:=Concatenation("gps",String(ind));
    ReadPackage( "primgrp", Concatenation( "data/", fname, ".g" ) );
  fi;
end);

BindGlobal("PrimGrpArtifactFilename",function(deg,nr)
  local filename;
  if deg <= 4095 then
    Error("This method is only for primitive groups of degree greater than 4095!");
  fi;
  filename:=Concatenation("PrimitiveGroups_", String(deg),"_", String(nr), ".g.gz");
  filename:=Filename(DirectoriesPackageLibrary("primgrp", "data/ExtendedPrimitiveGroupsData"), filename);
  return filename;
end);

BindGlobal("PRIMGrp",function(deg,nr)
  local filename,strm,r,l;
  if nr>PRIMLENGTHS[deg] then
    Error("There are only ",PRIMLENGTHS[deg]," groups of degree ",deg,"\n");
  fi;
  if deg > 4095 then
    filename:=PrimGrpArtifactFilename(deg,nr);
    if filename = fail then
      Error("Primitive group of degree ", deg, " with id ", nr, " not found! Note that primitive groups of degree 4096 to 8191 must be downloaded separately. They can be obtained from https://doi.org/10.5281/zenodo.10411366");
    fi;
    strm:=InputTextFile(filename);;
    r:=EvalString(ReadAll(strm));;
    CloseStream(strm);;
    if not "name" in RecNames(r) then
      r.name:="";
    fi;
    l:=[r.id, r.size, r.SimpleSolvable, r.ONanScottType, r.suborbits, r.transitivity, r.name, r.SocleType, r.generators];
    return l;
  fi;
  PrimGrpLoad(deg);
  l:=PRIMGRP[deg][nr];
  if IsFunction(l) then
    # An entry may be "lazy", that is, encoded in a function. We call such a
    # function the degree and index as arguments to produce the actual entry.
    # To avoid recomputing it, we store the the computed entry into `PRIMGRP`,
    # overwriting the function that produced it.
    l:=l(deg,nr);
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
  if deg <= 4095 then
    return true;
  elif deg <= 8191 then
    if PrimGrpArtifactFilename(deg,1) <> fail then
      return true;
    else
      Info(InfoWarning,1,"Note that primitive groups of degree 4096 to 8191 must be downloaded separately. They can be obtained from https://doi.org/10.5281/zenodo.10411366");
      return false;
    fi;
  else
    return false;
  fi;
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
  elif Length(l[9]) = 2 and l[9][1] = "4c" then
    # product action: the socle width in field 8 gives k, and the degree
    # its k-th root gives m
    k:=l[8][3];
    g:= PGProductAction4c(RootInt(deg,k), k, l[9][2]);
  elif l[4] = "1" and deg <= 4095 then
    # affine type groups described by matrices
    if Length(l[9]) > 0 then
      fac:= Factors(deg);
      mats:=List(l[9],i->ImmutableMatrix(GF(fac[1]),i));
      v:=Elements(GF(fac[1])^Length(fac));
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
