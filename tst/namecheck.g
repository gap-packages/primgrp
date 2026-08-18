#############################################################################
##
##  tst/namecheck.g
##
##  Checking the names in the library against StructureDescription.
##
##  Separate from tst/testutils.g because it binds constructors GAP does not
##  have -- AGL, ASL and friends -- which the names use and which would
##  otherwise have to be globals of the package.
##
##  Most names cannot be compared literally: a group has many correct
##  descriptions, and the library follows ATLAS conventions while
##  StructureDescription does not.  So the check is heuristic, and reports
##  ERROR only when a name evaluates to a group whose description differs,
##  WARNING when it cannot decide.  Some manual reading is expected.
##

if not IsBoundGlobal("AGL") then
  BindGlobal("AGL", {d,q} -> SemidirectProduct(GL(d, GF(q)), GF(q)^d));
fi;
if not IsBoundGlobal("ASL") then
  BindGlobal("ASL", {d,q} -> SemidirectProduct(SL(d, GF(q)), GF(q)^d));
fi;
if not IsBoundGlobal("L") then
  BindGlobal("L", {d,q} -> SimpleGroup("L", d, q));
fi;
if not IsBoundGlobal("ASigmaL") then
  ASigmaL := function(d, q)
    local G;
    G := SigmaL(d, q);
    return SemidirectProduct(G, FieldOfMatrixGroup(G)^DimensionOfMatrixGroup(G));
  end;
fi;
if not IsBoundGlobal("AGammaL") then
  AGammaL := function(d, q)
    local G;
    G := GammaL(d, q);
    return SemidirectProduct(G, FieldOfMatrixGroup(G)^DimensionOfMatrixGroup(G));
  end;
fi;

######################################################################
##
##  Evaluate a name and compare the result with <desc>.  Returns fail if the
##  name does not evaluate to a group.
##
ComparePrimGrpByStructure := function(name, desc)
  local n, g, r, brk, depth, c, i;
  n := ReplacedString( name, "A(",   "AlternatingGroup(" );
  n := ReplacedString( n,    "Alt(", "AlternatingGroup(" );
  n := ReplacedString( n,    "Aut(", "AutomorphismGroup(");
  n := ReplacedString( n,    "C(",   "CyclicGroup("      );
  n := ReplacedString( n,    "D(",   "DihedralGroup("    );
  n := ReplacedString( n,    "M(",   "MathieuGroup("     );
  n := ReplacedString( n,    "Q(",   "QuaternionGroup("  );
  n := ReplacedString( n,    "S(",   "SymmetricGroup("   );
  n := ReplacedString( n,    "Sym(", "SymmetricGroup("   );
  # Only evaluate a name that is exactly one constructor call, "Foo(...)".
  #
  # Anything else has to be refused rather than tried.  A malformed string
  # raises a *syntax* error, which CALL_WITH_CATCH does not catch: GAP drops
  # into a break loop and the run dies.  And a name ending ":2" is worse than
  # useless -- GAP reads the suffix as an option and discards it, so "M(22):2"
  # evaluates to M22 and gets reported as an order mismatch by a factor of
  # exactly 2, which is a bug in the check rather than in the data.
  if Length(n) = 0 or not IsAlphaChar(n[1]) or n[Length(n)] <> ')' then
    return fail;
  fi;
  depth := 0;
  for i in [1..Length(n)] do
    c := n[i];
    if c = '(' then
      depth := depth + 1;
    elif c = ')' then
      depth := depth - 1;
      if depth = 0 and i < Length(n) then
        return fail;                  # closed early, so more follows
      fi;
    elif depth = 0 and not (IsAlphaChar(c) or IsDigitChar(c)) then
      return fail;                    # anything but the leading identifier
    fi;
  od;
  if depth <> 0 then
    return fail;
  fi;
  # And the constructor has to exist.  PGammaU(3,5) is a perfectly well-formed
  # call to a function GAP does not have, and letting EvalString discover that
  # costs an error message and, with no terminal to fall back on, the run.
  i := Position(n, '(');
  if not IsBoundGlobal(n{[1..i-1]}) then
    return fail;
  fi;
  brk := BreakOnError;
  BreakOnError := false;
  r := CALL_WITH_CATCH( EvalString, [ n ] );
  BreakOnError := brk;
  if not r[1] or not IsGroup(r[2]) then
    return fail;                      # some names evaluate to a number
  fi;
  g := r[2];
  if StructureDescription( g ) = desc then
    return true;
  fi;
  # StructureDescription is not canonical, so two descriptions of the same
  # group can differ as text -- PGammaL(2,9) and PrimitiveGroup(10,7) are the
  # same group of order 1440 and describe differently.  Only a failed
  # isomorphism test makes it a real disagreement.
  return rec(differs := true, evaluated := g);
end;

##  "7:3" against "C7 : C3": drop the C from cyclic factors.
DropCfromCyclicGroups := function(str)
  local pos;
  str := ShallowCopy(str);
  repeat
    pos := Position(str, 'C');
    if pos <> fail and pos < Length(str) and IsDigitChar(str[pos+1]) then
      str := Concatenation(str{[1..pos-1]}, str{[pos+1..Length(str)]});
    elif pos <> fail then
      return str;
    fi;
  until pos = fail;
  return str;
end;

ComparePrimGrpByNames := function(name, desc)
  local n, d;
  n := ReplacedString(name, " ", "");
  n := ReplacedString(n, "Alt", "A");
  n := ReplacedString(n, "Sym", "S");
  d := ReplacedString(desc, " ", "");
  if (StartsWith(n,"A(") and EndsWith(n,")") and StartsWith(d,"A")) or
     (StartsWith(n,"C(") and EndsWith(n,")") and StartsWith(d,"C")) or
     (StartsWith(n,"M(") and EndsWith(n,")") and StartsWith(d,"M")) or
     (StartsWith(n,"S(") and EndsWith(n,")") and StartsWith(d,"S")) then
    n := ReplacedString(n, "(", "");
    n := ReplacedString(n, ")", "");
  elif StartsWith(n, "M_") then
    n := Concatenation("M", n{[3..Length(n)]});
  elif StartsWith(n,"D(") and EndsWith(n,")") and StartsWith(d,"D") then
    n := Concatenation("D", String(EvalString(n{[3..Length(n)-1]})));
  fi;
  return n = d or n = DropCfromCyclicGroups(d);
end;

##  Complaints about the names of one degree, as strings.
PrimGrpNameCheckDegree := function(deg)
  local out, i, g, desc, names, name, r;
  out := [];
  for i in [1..NrPrimitiveGroups(deg)] do
    g := PrimitiveGroup(deg, i);
    if not HasName(g) then
      continue;
    fi;
    desc := StructureDescription(g);
    names := SplitString(Name(g), "=");
    for name in names do
      if ComparePrimGrpByNames(name, desc) then
        continue;
      fi;
      r := ComparePrimGrpByStructure(name, desc);
      if r = true then
        continue;
      elif r = fail then
        Add(out, Concatenation("WARNING PrimitiveGroup(", String(deg), ",",
                 String(i), "): name ", name, " could not be evaluated; structure ", desc));
      elif Size(r.evaluated) <> Size(g) then
        Add(out, Concatenation("ERROR   PrimitiveGroup(", String(deg), ",",
                 String(i), "): name ", name, " has order ", String(Size(r.evaluated)),
                 ", group has ", String(Size(g))));
      elif IsomorphismGroups(r.evaluated, g) = fail then
        Add(out, Concatenation("ERROR   PrimitiveGroup(", String(deg), ",",
                 String(i), "): name ", name, " is not isomorphic to the group; structure ", desc));
      fi;
    od;
  od;
  return out;
end;
