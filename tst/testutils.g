######################################################################
#
# Checking the primitive groups of degree n
#
######################################################################
#
# 1. Consistency check: each group of degree n should act
# transitively and primitively on [1..n]
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

takeCfromcyclicgroup := function( str )
local pos;
repeat
  pos:=Position(str,'C');
  if pos <> fail then
    if IsDigitChar(str[pos+1]) then
      str := Concatenation(str{[1..pos-1]},str{[pos+1..Length(str)]});
    fi;
  fi;
until pos = fail;
return str;
end;

comparenames := function(name,desc)
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
return n=d or n=takeCfromcyclicgroup(d);
end;

findaffinegroups := function(type,name)
#
# For a finite field F
# AGL(d,F) = SemidirectProduct(GL(d, F), F^d)
# ASL(d,F) = SemidirectProduct(SL(d, F), F^d)
#
local n, old, new, pos, pos2, d, f, arguments;
n := ShallowCopy(name);
old := Concatenation( type, "(");
pos := PositionSublist( n, old );
if pos = 1 then
  pos2 := Position( n, ')' );
  arguments := n{[ Length(old)+1 .. pos2-1]};
  arguments := ReplacedString(arguments," ","");
  arguments := SplitString(arguments,",");
  new := Concatenation(
            "SemidirectProduct(",
            type{[2..3]}, # will be GL or SL
            "(", arguments[1], ",",
            "GF(", arguments[2], ")),",
            "GF(", arguments[2], ")^", arguments[1], ")"
            );
  n := Concatenation( new, n{[pos2+1 .. Length(n)]});
fi;
return n;
end;

findsimplegroups := function(type,name)
local n, old, new, pos;
n := ShallowCopy(name);
old := Concatenation( type, "(");
new := Concatenation( "SimpleGroup(\"",type,"\",");
pos := PositionSublist( name, old );
if pos = 1 then
    n := Concatenation(new, n{[Length(old)+1 .. Length(n)]});
fi;
return n;
end;

comparestructure := function (name, desc)
local n, g, r;
n := ReplacedString(name,"A(","AlternatingGroup(");
n := ReplacedString(n,"C(","CyclicGroup(");
n := ReplacedString(n,"M(","MathieuGroup(");
n := ReplacedString(n,"Q(","QuaternionGroup(");
n := ReplacedString(n,"S(","SymmetricGroup(");
n := findaffinegroups("AGL",n);
n := findaffinegroups("ASL",n);
n := findsimplegroups("L",n);
BreakOnError := false;
r := CALL_WITH_CATCH( EvalString, [ n ] );
BreakOnError := true;
if r[1] then
    g := r[2];
    if IsGroup(g) then
        return StructureDescription( g ) = desc;
    else
        return fail; # some may evaluate to a number
    fi;         
else
    return fail; 
fi;      
end;

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
                if comparenames( name, StructureDescription(g) ) then
                    continue;
                else 
                    r := comparestructure( name, StructureDescription(g) );
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

PrimGrpNamesCheck := function(i,j)
local k;
for k in [i..j] do
    PrimGrpNamesCheckDegree(k);
od;
Print("\n");
end;


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
