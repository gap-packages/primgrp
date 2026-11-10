#############################################################################
##
##  This file is part of the PrimGrp package for GAP.
##
##  Copyright of PrimGrp belongs to its developers; see the manual,
##  PackageInfo.g and the git log for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##  PrimitiveGroupsIterator, AllPrimitiveGroups and OnePrimitiveGroup.
##

##
#R  IsPrimGrpIterRep
##
DeclareRepresentation("IsPrimGrpIterRep",IsComponentObjectRep);

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
    Error("the arguments must be pairs of a function and a value");
  fi;
  for i in [1..l] do
    a:=arglis[2*i-1];
    b:=arglis[2*i];
    if not IsFunction(a) then
      Error("argument ",2*i-1," must be a function");
    fi;
    # degree and order values are used as integers below
    if a in [NrMovedPoints,Size,Order]
       and not (IsInt(b) or IsFunction(b)
                or (IsList(b) and (IsRange(b) or ForAll(b,IsInt)))) then
      Error("the value for ",NameFunction(a)," must be integers or a function");
    fi;
  od;
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

# needed by List, which iterates over a copy; the components are either
# immutable or replaced, never changed in place, so copying them suffices
InstallMethod(ShallowCopy,"primitive groups iterator",true,
  [IsPrimGrpIterRep and IsIterator],0,
function(it)
local new;
  new:=Objectify(NewType(IteratorsFamily,
                         IsIterator and IsPrimGrpIterRep and IsMutable),rec());
  new!.deg:=it!.deg;
  new!.degi:=it!.degi;
  new!.nr:=it!.nr;
  new!.prop:=it!.prop;
  new!.gut:=it!.gut;
  new!.next:=it!.next;
  return new;
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
