#############################################################################
##
##  dev/convert-ksets.g
##
##  Write the k-subset actions of Alt(n) and Sym(n) as calls.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-ksets.g
##
##  Everything such an entry records follows from n and k -- the order, the
##  Johnson suborbits, the transitivity, the socle -- so the entry becomes
##  PGAltOnSets(n,k) or PGSymOnSets(n,k) and the file keeps neither the
##  generators nor the rest.
##
##  Which of the two an entry is, is not decided by its name or its order.  For
##  each candidate the constructor is evaluated and compared against the entry;
##  and the group it names is built and measured as well, because field 9 stops
##  holding permutations and the group afterwards is built rather than read,
##  and so is a different conjugate.
##
##  Two fields are allowed to differ, and this is the point of saying so here
##  rather than in passing.
##
##  The name: three conventions are in the data at once -- "A(n)" up to degree
##  50, "Alt(n)" above it, and no name at all past degree 2555 -- and they
##  collapse to "Alt(n)" and "Sym(n)".  That is a change of data, and every one
##  is printed.
##
##  The order of the suborbits in field 5, which is compared as a set at its
##  one use in PrimitiveIdentification and so carries no information; the
##  entries hold it in whatever order their generators produced, and afterwards
##  it is sorted.
##
##  Nothing else about an entry is permitted to move.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

##  Which call, if any, reproduces this entry?  Returns the text to write, or
##  fail.
PRIMGRP_OnSetsCall := function(e, deg, report, resorted)
  local n, k, i, c, ctors, names;
  if not (IsList(e) and Length(e) = 9 and IsList(e[8]) and e[8][1] = "A"
          and IsInt(e[8][2]) and e[8][3] = 1) then
    return fail;
  fi;
  if IsFunction(e[9]) then
    return fail;                        # already a call
  fi;

  n := e[8][2];
  k := First([2..QuoInt(n,2)], i -> Binomial(n,i) = deg);
  if k = fail then
    return fail;                        # socle A_n, but not on its k-sets
  fi;

  ctors := [PGAltOnSets, PGSymOnSets];
  names := ["PGAltOnSets", "PGSymOnSets"];
  for i in [1,2] do
    c := ctors[i](n,k)(deg, e[1]);
    if c{[1,2,3,4,6]} = e{[1,2,3,4,6]} and Set(c[5]) = Set(e[5])
       and c[8] = e[8] then
      if c[7] <> e[7] then
        AddSet(report, [deg, e[1], e[7], c[7]]);
      fi;
      if c[5] <> e[5] then
        AddSet(resorted, [deg, e[1]]);
      fi;
      return rec(text := Concatenation(names[i], "(", String(n), ",",
                                       String(k), ")"),
                 inner := c[9][2], k := k);
    fi;
  od;
  return fail;
end;

##  Build the group the call names and measure it, rather than trusting that
##  the constructor's arithmetic and GAP's group agree.
PRIMGRP_CheckOnSets := function(r, e, deg)
  local g, subs;
  g := PGOnSetsGroup(r.inner, r.k);
  if NrMovedPoints(g) <> deg then
    return Concatenation("degree ", String(NrMovedPoints(g)));
  fi;
  if Size(g) <> e[2] then
    return Concatenation("order ", String(Size(g)));
  fi;
  if not IsPrimitive(g, [1..deg]) then
    return "not primitive";
  fi;
  if Transitivity(g, [1..deg]) <> e[6] then
    return Concatenation("transitivity ", String(Transitivity(g, [1..deg])));
  fi;
  subs := Collected(List(OrbitsDomain(Stabilizer(g, 1), [1..deg]), Length));
  subs := Filtered(subs, x -> x[1] <> 1);
  if SortedList(subs) <> SortedList(e[5]) then
    return Concatenation("suborbits ", String(subs));
  fi;
  return fail;
end;

##  A content line is one ending in a comma that begins with neither
##  "PRIMGRP[" nor "];".  That is what the normalisation pass bought.
PRIMGRP_ConvertFile := function(path, report, resorted)
  local lines, out, deg, line, e, r, why, n;
  lines := SplitString(StringFile(path), "\n");
  out := [];
  deg := fail;
  n := 0;
  for line in lines do
    if Length(line) = 0 then
      continue;
    fi;
    if PositionSublist(line, "PRIMGRP[") = 1 then
      deg := Int(line{[9..PositionSublist(line, "]") - 1]});
      Add(out, line);
    elif line{[1..2]} = "];" then
      Add(out, line);
    elif line[Length(line)] <> ',' then
      Add(out, line);
    elif PositionSublist(line, ",[\"A\",") = fail then
      Add(out, line);                   # cheap test before the dear one
    else
      e := EvalString(line{[1..Length(line)-1]});
      r := PRIMGRP_OnSetsCall(e, deg, report, resorted);
      if r = fail then
        Add(out, line);
      else
        why := PRIMGRP_CheckOnSets(r, e, deg);
        if why <> fail then
          Error("degree ", deg, " entry ", e[1], " would become ", r.text,
                " but that group has ", why);
        fi;
        Add(out, Concatenation(r.text, ","));
        n := n + 1;
      fi;
    fi;
  od;
  if n > 0 then
    FileString(path, Concatenation(JoinStringsWithSeparator(out, "\n"), "\n"));
  fi;
  return n;
end;

PRIMGRP_ConvertAll := function(dir)
  local files, f, n, report, resorted, r;
  files := Filtered(SortedList(DirectoryContents(dir)),
                    f -> Length(f) > 5 and f{[1..3]} = "gps"
                         and f{[Length(f)-1..Length(f)]} = ".g");
  n := 0;
  report := [];
  resorted := [];
  for f in files do
    n := n + PRIMGRP_ConvertFile(Concatenation(dir, "/", f), report, resorted);
  od;
  Print("CONVERTED ", n, " entries\n");
  Print("RESORTED ", Length(resorted), " suborbit lists\n");
  Print("RENAMED ", Length(report), " entries\n");
  for r in report do
    Print("NAME ", r[1], " ", r[2], " ", ViewString(r[3]), " -> ",
          ViewString(r[4]), "\n");
  od;
end;

PRIMGRP_ConvertAll(conv_dir);
QUIT;
