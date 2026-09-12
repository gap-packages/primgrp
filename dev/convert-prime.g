#############################################################################
##
##  dev/convert-prime.g
##
##  Write the affine entries of prime degree as calls.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-prime.g
##
##  For prime degree p the affine primitive groups are the subgroups of
##  AGL(1,p) containing the translations, one for each divisor d of p-1.
##  Everything the entry records follows from p and d, so it becomes ["Prime",d]
##  and the matrix goes away with the rest.
##
##  d is not guessed: it is the order divided by the degree, and the entry is
##  rewritten only if PGPrime(d) reproduces its first eight fields.  The group
##  is then built both ways and compared, because PGPrime chooses its own
##  generator of the subgroup of order d rather than keeping the stored one --
##  a different generator of the same cyclic group, which had better give the
##  same group.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

##  Build an affine group of prime degree the way the reader does.
PRIMGRP_PrimeGroup := function(deg, gens)
  local v, mats, perms, t;
  if Length(gens) = 0 then
    return Image(IsomorphismPermGroup(CyclicGroup(deg)));
  fi;
  v := AsSet(GF(deg)^1);
  mats := List(gens, m -> ImmutableMatrix(GF(deg), m));
  perms := List(mats, m -> Permutation(m, v, OnRight));
  t := First(v, x -> not IsZero(x));
  Add(perms, Permutation(t, v, function(a,b) return a+b; end));
  return Group(perms);
end;

##  The text for one entry, or fail.
PRIMGRP_ConvertPrime := function(e, deg)
  local d, c;
  if not (IsList(e) and Length(e) = 9 and e[4] = "1" and IsPrimeInt(deg)) then
    return fail;
  fi;
  # IsString([]) is true, so ask IsStringRep and a length: an entry whose
  # field 9 is a token is not an affine one, and one whose field 9 is the
  # empty list is -- the translations alone.
  if IsStringRep(e[9]) and Length(e[9]) > 0 then
    return fail;
  fi;

  d := e[2] / deg;
  if not IsInt(d) or (deg-1) mod d <> 0 then
    return fail;
  fi;

  c := PGPrime(d)(deg, e[1]);
  if c{[1..8]} <> e{[1..8]} then
    return fail;
  fi;

  if PRIMGRP_PrimeGroup(deg, e[9]) <> PRIMGRP_PrimeGroup(deg, c[9]) then
    Error("degree ", deg, " entry ", e[1], ": PGPrime(", d,
          ") builds a different group");
  fi;

  return Concatenation("[\"Prime\",", String(d), "]");
end;

##  A content line is one ending in a comma that begins with neither
##  "PRIMGRP[" nor "];".  That is what the normalisation pass bought.
PRIMGRP_ConvertFile := function(path)
  local lines, out, deg, line, e, text, n;
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
    elif not IsPrimeInt(deg) or PositionSublist(line, ",\"1\",") = fail then
      Add(out, line);                   # cheap test before the dear one
    else
      e := EvalString(line{[1..Length(line)-1]});
      text := PRIMGRP_ConvertPrime(e, deg);
      if text = fail then
        Add(out, line);
      else
        Add(out, Concatenation(text, ","));
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
