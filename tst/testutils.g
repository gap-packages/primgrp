######################################################################
#
# Checking the primitive groups of degree n
#
######################################################################
#
# 1. Consistency check: each group of degree n should act
# transitively and primitively on [1..n]
#

######################################################################
#
# Consistency check for a given degree
#
PrimGrpConsistencyCheckDegree := function(deg) 
    local i, g;
    for i in [1..NrPrimitiveGroups(deg)] do
        Print("Degree ", deg, " : ", i, "/", NrPrimitiveGroups(deg), "\r");
        g := PrimitiveGroup(deg,i);
        if MovedPoints(g) <> [1..deg] then
            Print("PrimitiveGroup(",deg,",",i,") is not acting on ", [1..deg], "\n");
        fi;
        if not IsTransitive(g,[1..deg]) then
            Print("PrimitiveGroup(",deg,",",i,") is not transitive on ", [1..deg], "\n");
        fi;
        if not IsPrimitive(g,[1..deg]) then
            Print("PrimitiveGroup(",deg,",",i,") is not primitive on ", [1..deg], "\n");
        fi;
    od;
    Print("                                                      \r");
end;

######################################################################
#
# Consistency check for a range of degrees
#
PrimGrpConsistencyCheck := function(i,j)
local k;
for k in [i..j] do
    PrimGrpConsistencyCheckDegree(k);
od;
Print("\n");
end;

######################################################################
#
# 2. Checking names of primitive groups
#

######################################################################
#
# Constructors for some groups
#
AGL := {d,q} -> SemidirectProduct(GL(d, GF(q)), GF(q)^d);
ASL := {d,q} -> SemidirectProduct(SL(d, GF(q)), GF(q)^d);
L   := {d,q} -> SimpleGroup("L", d, q);

ASigmaL := function(d, q)
    local G, V;
    G := SigmaL(d, q);
    V := FieldOfMatrixGroup(G) ^ DimensionOfMatrixGroup(G);
    return SemidirectProduct(G, V);
end;

AGammaL := function(d, q)
    local G, V;
    G := GammaL(d, q);
    V := FieldOfMatrixGroup(G) ^ DimensionOfMatrixGroup(G);
    return SemidirectProduct(G, V);
end;

if not IsBound(PSigmaL) then
  PSigmaL := function(d, q)
      local G, gens;
      if d = 1 then return TrivialGroup(IsPermGroup); fi;
      G := SigmaL(d, q);
      gens := GeneratorsOfGroup(G);
      if not IsPrimeInt(q) then
          gens := gens{[1..Length(gens)-1]}; # HACK: last generator is Frobenius
      fi;
      return G / Center(Subgroup(G, gens));
  end;
fi;

if not IsBound(PSigmaL) then
  PGammaL := function(d, q)
      local G, gens;
      if d = 1 then return TrivialGroup(IsPermGroup); fi;
      G := GammaL(d, q);
      gens := GeneratorsOfGroup(G);
      if not IsPrimeInt(q) then
          gens := gens{[1..Length(gens)-1]}; # HACK: last generator is Frobenius
      fi;
      return G / Center(Subgroup(G, gens));
  end;
fi;

######################################################################
#
# Try to evaluate the name of the primitive group and compare
# its structure description with the given one
#
ComparePrimGrpByStructure := function (name, desc)
local n, g, r;
n := ReplacedString( name, "A(",   "AlternatingGroup(" );
n := ReplacedString( n,    "Alt(", "AlternatingGroup(" );
n := ReplacedString( n,    "C(",   "CyclicGroup("      );
n := ReplacedString( n,    "D(",   "DihedralGroup("    );
n := ReplacedString( n,    "M(",   "MathieuGroup("     );
n := ReplacedString( n,    "Q(",   "QuaternionGroup("  );
n := ReplacedString( n,    "S(",   "SymmetricGroup("   );
n := ReplacedString( n,    "Sym(", "SymmetricGroup("   );
BreakOnError := false;
r := CALL_WITH_CATCH( EvalString, [ n ] );
BreakOnError := true;
if r[1] then
    g := r[2];
    if IsGroup(g) then # some names may evaluate to a number
        return StructureDescription( g ) = desc;
    else
        return fail;
    fi;
else
    return fail;
fi;
end;


######################################################################
#
# Remove C from cyclic groups in structure descriptions
#
# Deals with comparing names like "7:3" with description "C7 : C3"
#
DropCfromCyclicGroups := function( str )
local pos;
repeat
  pos:=Position(str,'C');
  if pos <> fail then
    # Is C followed by a digit?
    if IsDigitChar(str[pos+1]) then
      str := Concatenation(str{[1..pos-1]},str{[pos+1..Length(str)]});
    fi;
  fi;
until pos = fail;
return str;
end;

######################################################################
#
# Compare the name with the structure description
#
ComparePrimGrpByNames := function(name,desc)
local n, d;
n := ReplacedString(name," ","");
n := ReplacedString(n,"Alt","A");
n := ReplacedString(n,"Sym","S");
d := ReplacedString(desc," ","");
if StartsWith(n,"A(") and EndsWith(n,")") and StartsWith(d,"A") or
   StartsWith(n,"C(") and EndsWith(n,")") and StartsWith(d,"C") or
   StartsWith(n,"M(") and EndsWith(n,")") and StartsWith(d,"M") or
   StartsWith(n,"S(") and EndsWith(n,")") and StartsWith(d,"S") then
  n := ReplacedString(n,"(","");
  n := ReplacedString(n,")","");
elif StartsWith(n,"M_") then
  n := Concatenation( "M", n{[3..Length(n)]} );
elif StartsWith(n,"D(") and EndsWith(n,")") and StartsWith(d,"D") then
  n := Concatenation( "D", String(EvalString(n{[3..Length(n)-1]}) ) );
fi;
return n=d or n=DropCfromCyclicGroups(d);
end;


######################################################################
#
# Names check for a given degree
#
PrimGrpNamesCheckDegree := function(deg) 
    local i, g, r, names, name;
    for i in [1..NrPrimitiveGroups(deg)] do
        Print("Degree ", deg, " : ", i, "/", NrPrimitiveGroups(deg), "\r");
        g := PrimitiveGroup(deg,i);
        if HasName(g) then
            names := SplitString( Name(g), "=" );
            if Length(names) > 1 then
                 Print("WARNING : ", [deg,i], " has multiple names ", Name(g), "\n");
            fi;
            for name in names do
                if ComparePrimGrpByNames( name, StructureDescription(g) ) then
                    continue;
                else # if string transformations did not help
                    r := ComparePrimGrpByStructure( name, StructureDescription(g) );
                    if r = true then
                        continue;
                    elif r = false then
                        Print("ERROR : ", [deg,i], " : ", name, " is ", StructureDescription(g), "\n");
                    else       
                        Print("WARNING : ", [deg,i], " : ", name, " is ", StructureDescription(g), "\n");
                    fi;
                fi;        
            od;
        fi;
    od;
    Print("                                                      \r");
end;

######################################################################
#
# Names check for a range of degrees
#
PrimGrpNamesCheck := function(i,j)
local k;
for k in [i..j] do
    PrimGrpNamesCheckDegree(k);
od;
Print("\n");
end;

######################################################################
#
# Print all names of primitive groups
#
PrimGrpNames:=function()
local deg, i, g;
for deg in [1..4095] do
    for i in [1..NrPrimitiveGroups(deg)] do
        g := PrimitiveGroup(deg,i);
        if HasName(g) then
            Print( [deg, i ], " ", Name(g), "\n" );
        fi;
    od;
od;
end;
