#
# Checking the primitive groups of degree n
#
checkdegree := function(n) 
    local g;
    for g in AllPrimitiveGroups(DegreeOperation,n) do
        if MovedPoints(g) <> [1..n] or not IsTransitive(g,[1..n]) 
           or not IsPrimitive(g,[1..n]) then
            Error("Failure at ",g," degree ",n,"\n");
        fi;
    od;
end;
