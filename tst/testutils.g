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
