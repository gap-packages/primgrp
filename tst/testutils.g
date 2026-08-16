#
# Checking the primitive groups of degree n
#
checkdegree := function(n) 
    local g;
    for g in AllPrimitiveGroups(DegreeOperation,n) do
        if MovedPoints(g) <> [1..n] or not (IsNaturalAlternatingGroup(g) or IsNaturalSymmetricGroup(g) or IsPrimitive(g)) then
            Error("Failure at ",g," degree ",n,"\n");
        fi;
    od;
end;

#
# Check consistency of the cohorts data for the given degree
#
CheckCohortConsistency := function(deg)
    local tmp, cohort, i, type;
    tmp := Concatenation(List(COHORTS_PRIMITIVE_GROUPS[deg], x -> x[2]));
    Sort(tmp);
    if not IsSet(tmp) then
        Error("COHORTS_PRIMITIVE_GROUPS[",deg,"] lists the same group for multiple socle types");
    fi;
    if Length(tmp) < NrPrimitiveGroups(deg) then
        Error("COHORTS_PRIMITIVE_GROUPS[",deg,"] is leaving out some groups");
    elif Length(tmp) > NrPrimitiveGroups(deg) then
        Error("COHORTS_PRIMITIVE_GROUPS[",deg,"] covers too many (???) groups");
    fi;
    for cohort in COHORTS_PRIMITIVE_GROUPS[deg] do
        for i in cohort[2] do
            type := SocleTypePrimitiveGroup(PrimitiveGroup(deg,i));
            if type <> cohort[1] then
                Error("COHORTS_PRIMITIVE_GROUPS[",deg,"] mismatch with group ", i);
            fi;
        od;
    od;
end;

######################################################################
##
##  Stored versus recomputed.
##
##  PrimitiveGroup sets Size, Transitivity, ONanScottType, IsSimpleGroup,
##  IsSolvableGroup and SocleTypePrimitiveGroup from the data file rather than
##  computing them, and nothing has ever checked that the stored values are
##  right.  That is how the socle types of the natural A_n and S_n came to be
##  wrong for every degree from 2500 up (fixed in d80f6e7), and how 51 unitary
##  socles came to be labelled "^2A".
##
##  The recomputation has to happen on a *fresh* group, or it just reads the
##  attributes back.
##

BindGlobal("PrimGrpRecomputed", function(g)
  local h;
  h := Group(GeneratorsOfGroup(g));
  SetNrMovedPoints(h, NrMovedPoints(g));
  return h;
end);

# Collected suborbit lengths, as the data files store them.  Order varies --
# 4114 of the 24558 stored entries are not sorted -- so compare as sets, which
# is also what PrimitiveIdentification does.
BindGlobal("PrimGrpSuborbits", function(g, deg)
  return Set(Collected(List(Orbits(Stabilizer(g, 1), [2..deg]), Length)));
end);

##  A socle type as it appears in the alias list of
##  IsomorphismTypeInfoFiniteSimpleGroup: "L", [2,5] -> "L(2,5)".
BindGlobal("PrimGrpSocleToken", function(st)
  local par;
  if IsList(st.parameter) then
    par := JoinStringsWithSeparator(List(st.parameter, String), ",");
  else
    par := String(st.parameter);
  fi;
  return Concatenation(st.series, "(", par, ")");
end);

##  Do two socle types describe the same group?
##
##  Not just a record comparison: the simple groups have exceptional
##  isomorphisms, so PSL(2,5) is equally correctly "L", [2,5] and "A", 5, and
##  the library and IsomorphismTypeInfoFiniteSimpleGroup need not pick the same
##  one.  The `name' field of the computed type lists every alias --
##  "A(5) ~ A(1,4) = L(2,4) ~ ... ~ A(1,5) = L(2,5) ~ ..." -- so a stored type
##  agrees if its token appears there.
BindGlobal("PrimGrpSameSocleType", function(stored, computed)
  if stored.width <> computed.width then
    return false;
  fi;
  if stored.series = computed.series and stored.parameter = computed.parameter then
    return true;
  fi;
  return IsBound(computed.name) and
         PositionSublist(computed.name, PrimGrpSocleToken(stored)) <> fail;
end);

##  Check one group.  Returns a list of complaints, empty if all is well.
##  <checks> selects what to recompute; the socle is much the most expensive.
BindGlobal("PrimGrpCheckGroup", function(deg, nr, checks)
  local out, g, h, note, s, t;

  out := [];
  note := function(what, stored, found)
    Add(out, Concatenation("PrimitiveGroup(", String(deg), ",", String(nr),
             ") ", what, ": stored ", String(stored), ", computed ", String(found)));
  end;

  g := PrimitiveGroup(deg, nr);

  if MovedPoints(g) <> [1..deg] and not IsTrivial(g) then
    note("moved points", [1..deg], MovedPoints(g));
  fi;
  if not (IsNaturalAlternatingGroup(g) or IsNaturalSymmetricGroup(g)) then
    if not IsTransitive(g, [1..deg]) then
      note("transitive", true, false);
    elif not IsPrimitive(g, [1..deg]) then
      note("primitive", true, false);
    fi;
  fi;

  h := PrimGrpRecomputed(g);

  if "size" in checks and Size(h) <> Size(g) then
    note("Size", Size(g), Size(h));
  fi;
  if "transitivity" in checks then
    t := Transitivity(h, [1..deg]);
    if t <> Transitivity(g) then note("Transitivity", Transitivity(g), t); fi;
  fi;
  if "simple" in checks then
    if IsSimpleGroup(h) <> IsSimpleGroup(g) then
      note("IsSimpleGroup", IsSimpleGroup(g), IsSimpleGroup(h));
    fi;
    if IsSolvableGroup(h) <> IsSolvableGroup(g) then
      note("IsSolvableGroup", IsSolvableGroup(g), IsSolvableGroup(h));
    fi;
  fi;
  if "suborbits" in checks then
    s := PrimGrpSuborbits(g, deg);
    t := Set(PRIMGrp(deg, nr)[5]);
    if s <> t then note("suborbits", t, s); fi;
  fi;
  if "onanscott" in checks and ONanScottType(h) <> ONanScottType(g) then
    note("ONanScottType", ONanScottType(g), ONanScottType(h));
  fi;
  if "socle" in checks then
    s := SocleTypePrimitiveGroup(g);
    t := SocleTypePrimitiveGroup(h);
    if not PrimGrpSameSocleType(s, t) then
      note("SocleTypePrimitiveGroup",
           [s.series, s.parameter, s.width], [t.series, t.parameter, t.width]);
    fi;
  fi;

  return out;
end);

BindGlobal("PrimGrpCheckDegree", function(deg, checks)
  return Concatenation(List([1..NrPrimitiveGroups(deg)],
                            nr -> PrimGrpCheckGroup(deg, nr, checks)));
end);

BindGlobal("PrimGrpCheckAll", [ "size", "transitivity", "simple",
                                "suborbits", "onanscott", "socle" ]);
