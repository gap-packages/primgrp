#############################################################################
##
##  This file is part of the PrimGrp package for GAP.
##
##  Copyright of PrimGrp belongs to its developers; see the manual,
##  PackageInfo.g and the git log for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##  Loading the library entries, and building a group from one.
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

# Refuse a degree outside the library, naming the function the user called.
BindGlobal("PRIMGRP_CheckDegree", function(fname, deg)
  if not IsPosInt(deg) then
    Error(fname, ": <deg> must be a positive integer");
  fi;
  if deg > Maximum(PRIMRANGE) then
    Error(fname, ": <deg> must be at most ", Maximum(PRIMRANGE));
  fi;
end);

InstallGlobalFunction(NrPrimitiveGroups, function(deg)
  PRIMGRP_CheckDegree("NrPrimitiveGroups", deg);
  if not IsBound(PRIMLENGTHS[deg]) then
    PrimGrpLoad(deg);
  fi;
  return PRIMLENGTHS[deg];
end);

InstallGlobalFunction(PrimitiveGroupsAvailable,function(deg)
  return deg in PRIMRANGE;
end);

InstallGlobalFunction( PrimitiveGroup, function(deg,num)
local l,g,fac,mats,perms,v,t,filename,strm,r,dim,q,k,nr;

  PRIMGRP_CheckDegree("PrimitiveGroup", deg);
  if not IsPosInt(num) then
    Error("PrimitiveGroup: <nr> must be a positive integer");
  fi;
  nr:=NrPrimitiveGroups(deg);
  if num > nr then
    Error("PrimitiveGroup: there are only ",nr,
          " primitive groups of degree ",deg);
  fi;

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
    g:= PGPslExtendedGroup(l[9]);
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
