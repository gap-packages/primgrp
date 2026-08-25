#############################################################################
##
##  dev/normalise-data.g
##
##  Rewrite data/gps*.g so that every entry is one line.
##
##  The files currently wrap an entry across lines and indent it into the
##  degree's list, so every tool that edits them has had to walk brackets over
##  the whole file.  Six such tools have been written for this work and three
##  of them had bugs that were purely about the layout -- one found 87 of 406
##  entries and silently skipped the rest.  After this, a content line is one
##  that ends in a comma and begins with neither "PRIMGRP[" nor "];".
##
##  Nothing about the data changes.  This reads the files with GAP's own
##  reader and prints the values back, so the output cannot disagree with the
##  input about what an entry is -- which is the reason to do it in GAP rather
##  than by editing text.
##
##  Run it as
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'norm_dir:="data";;' \
##          dev/normalise-data.g
##
##  which rewrites every data/gps*.g in place and checks each one as it goes.
##  Nothing has to be backed up first: the entries are kept in memory, the file
##  is written, and then read again and compared against what was kept.  An
##  in-place rewrite that needed a pristine copy beside it would be a step
##  someone has to remember, and this has none.
##
##  It has to happen now.  An entry may already be a function, and PGAlt is a
##  global whose name NameFunction can recover; once entries are calls like
##  PGPsl(2,5) they evaluate to anonymous closures and no round trip can write
##  them back.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);

##  p if <o> is a nested list bottoming out in elements of the prime field
##  GF(p), fail otherwise.  Matrices are stored that way and are most of the
##  bulk: "0*Z(2)" and "Z(2)^0" become "0" and "1" under one Z(p)^0 factor.
PRIMGRP_FFEChar := function(o)
  local c, d, x;
  if IsFFE(o) then
    if DegreeFFE(o) <> 1 then
      return fail;                  # IntFFE is not the inverse off the
    fi;                             # prime field
    return Characteristic(o);
  fi;
  if not IsList(o) or IsEmpty(o) or IsString(o) then
    return fail;
  fi;
  c := fail;
  for x in o do
    d := PRIMGRP_FFEChar(x);
    if d = fail then
      return fail;
    elif c = fail then
      c := d;
    elif c <> d then
      return fail;
    fi;
  od;
  return c;
end;

PRIMGRP_IntStructure := function(o)
  if IsFFE(o) then
    return String(IntFFE(o));
  fi;
  return Concatenation("[", JoinStringsWithSeparator(
           List(o, PRIMGRP_IntStructure), ","), "]");
end;

PRIMGRP_Compact := function(o)
  local p;
  if IsFunction(o) then
    return NameFunction(o);
  fi;
  if IsStringRep(o) then
    # ViewString quotes and escapes; the name of PrimitiveGroup(625,657)
    # contains a backslash, and writing it unescaped is exactly the fault
    # that made it read back as "4 <chr 3>irc"
    return ViewString(o);
  fi;
  p := PRIMGRP_FFEChar(o);
  if p <> fail then
    return Concatenation("Z(", String(p), ")^0*", PRIMGRP_IntStructure(o));
  fi;
  if IsList(o) then
    return Concatenation("[", JoinStringsWithSeparator(
             List(o, PRIMGRP_Compact), ","), "]");
  fi;
  return String(o);
end;

PRIMGRP_Clear := function()
  local i;
  for i in [1..Length(PRIMGRP)] do
    Unbind(PRIMGRP[i]);
  od;
end;

PRIMGRP_Degrees := function()
  return Filtered([1..Length(PRIMGRP)], i -> IsBound(PRIMGRP[i]));
end;

##  Rewrite one file and check it: the entries are held while the file is
##  replaced, then read back and compared as values.  Not as printed text --
##  printing both sides would hide a fault in the printer by making it twice.
PRIMGRP_NormaliseFile := function(path)
  local degs, keep, d, out, e, n, bad, j, after;
  PRIMGRP_Clear();
  Read(path);                       # so the degrees now bound are this file's
  degs := PRIMGRP_Degrees();
  keep := [];
  for d in degs do
    keep[d] := PRIMGRP[d];
  od;

  n := 0;
  out := OutputTextFile(path, false);
  SetPrintFormattingStatus(out, false);
  for d in degs do
    PrintTo(out, "PRIMGRP[", d, "]:=[\n");
    for e in keep[d] do
      PrintTo(out, PRIMGRP_Compact(e), ",\n");
      n := n + 1;
    od;
    PrintTo(out, "];\n");
  od;
  CloseStream(out);

  PRIMGRP_Clear();
  Read(path);
  after := PRIMGRP_Degrees();
  bad := 0;
  if after <> degs then
    Print("BAD ", path, " defines different degrees after rewriting\n");
    return fail;
  fi;
  for d in degs do
    if Length(keep[d]) <> Length(PRIMGRP[d]) then
      Print("BAD ", path, " degree ", d, " changed length\n");
      bad := bad + 1;
    else
      for j in [1..Length(keep[d])] do
        if keep[d][j] <> PRIMGRP[d][j] then
          bad := bad + 1;
          if bad <= 3 then
            Print("BAD ", path, " degree ", d, " entry ", j, " differs\n");
          fi;
        fi;
      od;
    fi;
  od;
  return rec(degrees := Length(degs), entries := n, bad := bad);
end;

PRIMGRP_NormaliseAll := function(dir)
  local files, f, r, degs, n, bad;
  files := Filtered(SortedList(DirectoryContents(dir)),
                    f -> Length(f) > 5 and f{[1..3]} = "gps"
                         and f{[Length(f)-1..Length(f)]} = ".g");
  degs := 0;
  n := 0;
  bad := 0;
  for f in files do
    r := PRIMGRP_NormaliseFile(Concatenation(dir, "/", f));
    if r = fail then
      Error("normalising ", f, " changed which degrees it defines");
    fi;
    degs := degs + r.degrees;
    n := n + r.entries;
    bad := bad + r.bad;
  od;
  Print("NORMALISED ", Length(files), " files, ", degs, " degrees, ",
        n, " entries; ", bad, " differ after reading back\n");
  if bad <> 0 then
    Error(bad, " entries changed value, which this must never do");
  fi;
end;

PRIMGRP_NormaliseAll(norm_dir);
QUIT;
