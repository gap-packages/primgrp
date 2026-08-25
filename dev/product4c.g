#############################################################################
##
##  dev/product4c.g
##
##  Write the type 4c entries as elements of Sym(m) wreath Sym(k).
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/product4c.g
##
##  A type 4c group acts on the m^k tuples over a set of size m.  Its socle is
##  T^k, and the orbits of all the factors but the i-th are the fibres of the
##  i-th coordinate; that recovers the identification of the points with
##  tuples, after which each generator reads off as (p_1,...,p_k; sigma).
##  Storing those costs k permutations of degree m and one of degree k in
##  place of a permutation of degree m^k.
##
##  The relabelling of the points is a bijection, so the group built back from
##  the tuples is the original conjugated by it -- conjugate, not equal.  The
##  check below is on that conjugate, and so is an identity rather than a
##  comparison of invariants; and because a conjugate has the same order,
##  transitivity, suborbits and socle, fields 1 to 8 keep their meaning
##  without being re-measured.
##
##  It is made generator by generator rather than between the two groups.
##  Each generator is decomposed on its own, so that is the statement actually
##  being proved, and it is the cheaper one: comparing two groups of degree
##  2116 and order 10^115 costs minutes where comparing their generators costs
##  nothing.
##
##  Splitting the socle into its factors is the expensive step.  Asking for
##  MinimalNormalSubgroups is hopeless on the large ones -- degree 1089 entry 6,
##  whose socle is A(33)^2, does not finish in twenty minutes -- so the factors
##  are found by chance instead, in four seconds there.  See PRIMGRP_Components
##  below; nothing rests on the guess being right, because a wrong factor gives
##  the wrong number of fibres, and past that the rebuilt group fails to be the
##  relabelled original.
##
##  Even so a few entries may be worse than the rest, so this is built to be
##  killed:
##
##    * it announces each entry before attempting it, so whatever a kill
##      interrupts is named in the log,
##    * it rewrites the file after every conversion, so a kill costs one entry
##      rather than a file, and
##    * an entry already converted is recognised and left alone, so a re-run
##      resumes.
##
##  Pass conv_skip, a list of [degree, number], to leave the known-expensive
##  ones for later.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/product4c.g");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

if not IsBound(conv_skip) then
  conv_skip := [];
fi;

##  Print without the whitespace that String puts in lists.
PRIMGRP_Compact := function(o)
  if IsStringRep(o) then
    return ViewString(o);
  fi;
  if IsList(o) then
    return Concatenation("[", JoinStringsWithSeparator(
             List(o, PRIMGRP_Compact), ","), "]");
  fi;
  return String(o);
end;

##  The text for one entry, or fail if it is not a type 4c one to convert.
PRIMGRP_Convert4c := function(e, deg)
  local k, G, r;
  if not (IsList(e) and Length(e) = 9 and e[4] = "4c") then
    return fail;
  fi;
  if Length(e[9]) = 2 and e[9][1] = "4c" then
    return fail;                        # already converted
  fi;
  if [deg, e[1]] in conv_skip then
    Print("LEFT ", deg, " ", e[1], "\n");
    return fail;
  fi;

  Print("TRY ", deg, " ", e[1], "\n");
  k := e[8][3];
  G := GroupByGenerators(e[9], ());
  r := PRIMGRP_DecomposeProductAction(G, k);
  if IsBound(r.err) then
    Print("SKIP ", deg, " ", e[1], " ", r.err, "\n");
    return fail;
  fi;
  if not PRIMGRP_CheckProductAction(G, r) then
    Error("degree ", deg, " entry ", e[1],
          ": a generator does not rebuild as the relabelled original");
  fi;

  e := ShallowCopy(e);
  e[9] := ["4c", r.els];
  return PRIMGRP_Compact(e);
end;

##  A content line is one ending in a comma that begins with neither
##  "PRIMGRP[" nor "];".  That is what the normalisation pass bought.
PRIMGRP_ConvertFile := function(path)
  local lines, out, deg, i, line, e, text, n, t, flush;

  lines := SplitString(StringFile(path), "\n");
  out := [];
  deg := fail;
  n := 0;

  # Write out what has been converted so far followed by the untouched tail,
  # which is exactly what a resumed run expects to find.  It goes through a
  # temporary file, so a kill during the write cannot leave a half-written
  # data file behind.
  flush := function(i)
    local all;
    all := Concatenation(out, lines{[i+1..Length(lines)]});
    FileString(Concatenation(path, ".new"),
               Concatenation(JoinStringsWithSeparator(all, "\n"), "\n"));
    Exec(Concatenation("mv ", path, ".new ", path));
  end;

  for i in [1..Length(lines)] do
    line := lines[i];
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
    elif PositionSublist(line, ",\"4c\",") = fail then
      Add(out, line);                   # cheap test before the dear one
    else
      e := EvalString(line{[1..Length(line)-1]});
      t := Runtime();
      text := PRIMGRP_Convert4c(e, deg);
      if text = fail then
        Add(out, line);
      else
        Print("DONE ", deg, " ", e[1], " ", Runtime()-t, " ms\n");
        Add(out, Concatenation(text, ","));
        n := n + 1;
        flush(i);
      fi;
    fi;
  od;
  return n;
end;

PRIMGRP_ConvertAll := function(dir)
  local files, f, n;
  files := Filtered(SortedList(DirectoryContents(dir)),
                    f -> Length(f) > 5 and f{[1..3]} = "gps"
                         and f{[Length(f)-1..Length(f)]} = ".g");
  n := 0;
  for f in files do
    n := n + PRIMGRP_ConvertFile(Concatenation(dir, "/", f));
  od;
  Print("CONVERTED ", n, " entries\n");
end;

PRIMGRP_ConvertAll(conv_dir);
QUIT;
