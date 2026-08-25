#############################################################################
##
##  dev/convert-projective.g
##
##  Write the natural projective actions of the L series as calls.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-projective.g
##
##  PSL, PGL, PSigmaL and PGammaL on the points of PG(dim-1,q).  Everything
##  the entry records follows from dim and q, so the entry becomes the call.
##
##  Which of the four an entry is, is not decided by its name.  For each
##  candidate the constructor is evaluated and its first eight fields compared
##  against the entry's: an entry is rewritten only as something already known
##  to reproduce it.  Field 9 is the exception and is the point -- it holds
##  explicit permutations for most of these, and afterwards names the group
##  instead, so the group is built rather than read and is a different
##  conjugate.  That is the one thing here that is a change and not a
##  rewriting, and it is why each converted group is also constructed and
##  measured below.
##
##  Entries the four do not reproduce are left alone and counted.  Fourteen
##  are expected: intermediate subgroups between PSL and PGammaL which none of
##  the four names.
##
##  The file is edited as lines, since one entry is one line, but a candidate
##  line is turned back into its value with EvalString before anything is
##  decided about it.  The cheap text test first -- does the line mention an
##  L socle -- keeps that off the 15000 lines that could never match, some of
##  which are fifty thousand characters long.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

PRIMGRP_ProjectiveNames := ["psl", "pgl", "psigmal", "pgammal"];

##  the four constructors, in the same order
PRIMGRP_ProjectiveCtors := [PGPsl, PGPgl, PGPsigmaL, PGPgammaL];

##  Which call, if any, reproduces this entry?  Returns the text to write, or
##  fail.
PRIMGRP_ProjectiveCall := function(e, deg)
  local dim, q, i, c;
  if not IsList(e) or Length(e) <> 9 then
    return fail;
  fi;
  if not (IsList(e[8]) and Length(e[8]) = 3 and e[8][1] = "L"
          and IsList(e[8][2]) and Length(e[8][2]) = 2) then
    return fail;
  fi;
  dim := e[8][2][1];
  q := e[8][2][2];
  if (q^dim-1)/(q-1) <> deg then
    return fail;                    # socle L(dim,q) but not on its points
  fi;
  for i in [1..4] do
    c := PRIMGRP_ProjectiveCtors[i](dim, q)(deg, e[1]);
    if c{[1..8]} = e{[1..8]} then
      return rec(text := Concatenation(
                   ["PGPsl", "PGPgl", "PGPsigmaL", "PGPgammaL"][i],
                   "(", String(dim), ",", String(q), ")"),
                 dim := dim, q := q, which := PRIMGRP_ProjectiveNames[i]);
    fi;
  od;
  return fail;
end;

##  Build the group the converted entry names and measure it, rather than
##  trusting that the constructor's arithmetic and GAP's group agree.
PRIMGRP_CheckProjective := function(r, e, deg)
  local g, s, subs;
  if r.which = "psl" then
    g := PSL(r.dim, r.q);
  elif r.which = "pgl" then
    g := PGL(r.dim, r.q);
  elif r.which = "psigmal" then
    g := PSigmaL(r.dim, r.q);
  else
    g := PGammaL(r.dim, r.q);
  fi;
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
PRIMGRP_ConvertFile := function(path, report)
  local lines, out, deg, line, body, e, r, why, n, left, m;
  lines := SplitString(StringFile(path), "\n");
  out := [];
  deg := fail;
  n := 0;
  left := 0;
  for line in lines do
    if Length(line) = 0 then
      continue;
    fi;
    m := PositionSublist(line, "PRIMGRP[");
    if m = 1 then
      deg := Int(line{[9..PositionSublist(line, "]") - 1]});
      Add(out, line);
    elif line{[1..2]} = "];" then
      Add(out, line);
    elif line[Length(line)] <> ',' then
      Add(out, line);
    elif PositionSublist(line, ",[\"L\",[") = fail then
      Add(out, line);               # cheap test before the dear one
    else
      body := line{[1..Length(line)-1]};
      e := EvalString(body);
      r := PRIMGRP_ProjectiveCall(e, deg);
      if r = fail then
        Add(out, line);
        left := left + 1;
      else
        why := PRIMGRP_CheckProjective(r, e, deg);
        if why <> fail then
          Error("degree ", deg, " entry ", e[1], " would become ", r.text,
                " but that group has ", why);
        fi;
        Add(out, Concatenation(r.text, ","));
        n := n + 1;
        AddSet(report, [deg, e[1], r.which]);
      fi;
    fi;
  od;
  if n > 0 then
    FileString(path, JoinStringsWithSeparator(out, "\n"));
    AppendTo(path, "\n");
  fi;
  return [n, left];
end;

PRIMGRP_ConvertAll := function(dir)
  local files, f, r, n, left, report;
  files := Filtered(SortedList(DirectoryContents(dir)),
                    f -> Length(f) > 5 and f{[1..3]} = "gps"
                         and f{[Length(f)-1..Length(f)]} = ".g");
  n := 0;
  left := 0;
  report := [];
  for f in files do
    r := PRIMGRP_ConvertFile(Concatenation(dir, "/", f), report);
    n := n + r[1];
    left := left + r[2];
  od;
  Print("CONVERTED ", n, " entries; ", left,
        " with an L socle on its points that none of the four reproduces\n");
  Print("BY KIND ", Collected(List(report, x -> x[3])), "\n");
end;

PRIMGRP_ConvertAll(conv_dir);
QUIT;
