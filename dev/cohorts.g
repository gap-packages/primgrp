#############################################################################
##
##  dev/cohorts.g
##
##  Recompute COHORTS_PRIMITIVE_GROUPS from the data files and compare it with
##  the stored array.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" dev/cohorts.g
##      gap -q -b -A --quitonbreak -l "ROOT;" \
##          -c 'cohorts_degrees:=[2..8191];; cohorts_out:="cohorts.out";;' \
##          dev/cohorts.g
##
##  cohorts_degrees defaults to the degrees the stored array covers.  If
##  cohorts_out is set, the recomputed array is written there in the notation
##  of lib/cohorts.grp, one cohort per line.
##
##  A cohort is the set of groups of one degree that share a socle type.  The
##  socle type is field 8 of a PRIMGRP entry, so no group is built.  An entry
##  stored as a description is evaluated by PRIMGrp; for a projective group on
##  k-sets that computes suborbits in a small permutation group.
##
##  Cohorts are listed in the order of their first group, as in the stored
##  array.  A degree whose cohorts agree with the stored ones only up to that
##  order is reported separately.
##
##  Dixon and Mortimer group by the socle as a permutation group, not by its
##  isomorphism type: inequivalent actions of one socle are separate cohorts,
##  as PSU(3,3) on 63 points in [1,2] and [3,4].  Socle types cannot separate
##  those, so a stored array that splits a recomputed cohort, and agrees with
##  it otherwise, is reported as a split, not as a difference.
##

LoadPackage("primgrp");
SetPrintFormattingStatus("*stdout*",false);

if not IsBound(cohorts_degrees) then
  cohorts_degrees:=[2..Length(COHORTS_PRIMITIVE_GROUPS)];
fi;

# The cohorts of degree deg, as in COHORTS_PRIMITIVE_GROUPS[deg].
CohortsOfDegree:=function(deg)
  local types,cohorts,nr,s,type,pos;
  types:=[];
  cohorts:=[];
  for nr in [1..NrPrimitiveGroups(deg)] do
    s:=PRIMGrp(deg,nr)[8];
    type:=rec(series:=s[1], parameter:=s[2], width:=s[3]);
    pos:=Position(types,type);
    if pos = fail then
      Add(types,type);
      Add(cohorts,[type,[nr]]);
    else
      Add(cohorts[pos][2],nr);
    fi;
  od;
  return cohorts;
end;

# Order-independent form of a cohort list.
CohortsNormalised:=cohorts->Set(cohorts,c->[c[1],Set(c[2])]);

# The cohorts with the stored ones of one socle type merged.
CohortsMerged:=function(cohorts)
  local types;
  types:=Set(cohorts,c->c[1]);
  return Set(types,t->[t,Union(List(Filtered(cohorts,c->c[1]=t),c->c[2]))]);
end;

# A cohort in the notation of lib/cohorts.grp.
CohortLine:=function(c)
  local par;
  par:=c[1].parameter;
  if IsString(par) then
    par:=Concatenation("\"",par,"\"");
  else
    par:=String(par);
  fi;
  return Concatenation("[ rec( series:= \"",c[1].series,"\", parameter:= ",
    par,", width:= ",String(c[1].width),"),",String(c[2]),"]");
end;

CohortString:=c->Concatenation(c[1].series," ",String(c[1].parameter),
  " width ",String(c[1].width),": ",String(c[2]));

cohorts_new:=[];
cohorts_differ:=[];
cohorts_reordered:=[];
cohorts_split:=[];
cohorts_unstored:=[];

for deg in cohorts_degrees do
  new:=CohortsOfDegree(deg);
  cohorts_new[deg]:=new;

  # PRIMGrp keeps every entry it loads; the next degree needs none of them.
  Unbind(PRIMGRP[deg]);

  if not IsBound(COHORTS_PRIMITIVE_GROUPS[deg]) then
    Add(cohorts_unstored,deg);
    continue;
  fi;
  old:=COHORTS_PRIMITIVE_GROUPS[deg];
  if old = new then
    continue;
  fi;
  if CohortsNormalised(old) = CohortsNormalised(new) then
    Add(cohorts_reordered,deg);
    continue;
  fi;

  if CohortsMerged(old) = CohortsNormalised(new) then
    Add(cohorts_split,deg);
    Print("split ",deg,":");
    for c in Difference(CohortsNormalised(old),CohortsNormalised(new)) do
      Print("  ",CohortString(c),";");
    od;
    Print("\n");
    continue;
  fi;

  Add(cohorts_differ,deg);
  Print("degree ",deg,":\n");
  for c in Difference(CohortsNormalised(old),CohortsNormalised(new)) do
    Print("  stored:     ",CohortString(c),"\n");
  od;
  for c in Difference(CohortsNormalised(new),CohortsNormalised(old)) do
    Print("  recomputed: ",CohortString(c),"\n");
  od;
od;

Print("degrees checked: ",Length(cohorts_degrees),"\n");
Print("degrees differing: ",Length(cohorts_differ)," ",cohorts_differ,"\n");
Print("degrees differing only in cohort order: ",Length(cohorts_reordered),
      " ",cohorts_reordered,"\n");
Print("degrees with a socle type split into several cohorts: ",
      Length(cohorts_split)," ",cohorts_split,"\n");
Print("degrees not stored: ",Length(cohorts_unstored),"\n");

if IsBound(cohorts_out) then
  out:=OutputTextFile(cohorts_out,false);
  SetPrintFormattingStatus(out,false);
  PrintTo(out,"COHORTS_PRIMITIVE_GROUPS:= [];\n");
  for deg in cohorts_degrees do
    PrintTo(out,"COHORTS_PRIMITIVE_GROUPS[",deg,"]:= [\n",
      JoinStringsWithSeparator(List(cohorts_new[deg],CohortLine),",\n"),
      "];\n");
  od;
  CloseStream(out);
  Print("written to ",cohorts_out,"\n");
fi;

QUIT_GAP(0);
