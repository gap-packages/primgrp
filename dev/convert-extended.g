#############################################################################
##
##  dev/convert-extended.g
##
##  Convert the degree 4096-8191 data from one gzipped record per group into
##  ordinary data/gps*.g files.
##
##  The point of the exercise is the affine groups: 19571 of the 29776 groups
##  are O'Nan-Scott type 1 and are stored as permutations of degree up to 8191,
##  which is 94.7% of the 1 GB.  primgrp already has a matrix encoding for
##  exactly those (lib/primitiv.gi:121-134); the only reason it was not used is
##  the `deg <= 4095' guard.  The point labelling matches, so the conversion is
##  lossless -- PRIMGRP_CheckAffineRoundTrip insists on it group by group.
##
##  Usage:
##    gap -q -b dev/convert-extended.g   (set PRIMGRP_ExtDir/PRIMGRP_OutDir first)
##

if not IsBound(PRIMGRP_ExtDir) then
  PRIMGRP_ExtDir := "ExtendedPrimitiveGroupsData";
fi;
if not IsBound(PRIMGRP_OutDir) then
  PRIMGRP_OutDir := "converted";
fi;
if not IsBound(PRIMGRP_Degrees) then
  PRIMGRP_Degrees := [4096..8191];
fi;

LoadPackage("primgrp");
Read("dev/affine.g");
Read("dev/generic.g");
Read("dev/serialize.g");

# Read one legacy per-group file.  This is the last use of EvalString; the
# whole point is to stop shipping data that needs it (issue #60).
BindGlobal("PRIMGRP_ReadLegacy", function(deg, nr)
  local name, strm, r;
  name := Concatenation(PRIMGRP_ExtDir, "/PrimitiveGroups_",
                        String(deg), "_", String(nr), ".g.gz");
  strm := InputTextFile(name);
  if strm = fail then
    Error("cannot read ", name);
  fi;
  r := EvalString(ReadAll(strm));
  CloseStream(strm);
  return r;
end);

# Legacy record -> the 9-element PRIMGRP entry.
BindGlobal("PRIMGRP_LegacyEntry", function(r, deg, nr)
  local size, name;
  size := r.size;
  if IsString(size) then           # "Factorial(8191)", "Factorial(5222)/2"
    size := EvalString(size);      # stored as a *string*, so lib/primitiv.gi:202
  fi;                              # never matches it against an integer Size
  name := "";
  if IsBound(r.name) then
    name := r.name;
  fi;
  return [ r.id, size, r.SimpleSolvable, r.ONanScottType, r.suborbits,
           r.transitivity, name, r.SocleType, r.generators ];
end);

BindGlobal("PRIMGRP_ConvertDegree", function(deg)
  local n, entries, nr, r, e, gen, diff, G, rec_, mats, enum, stats;
  n := NrPrimitiveGroups(deg);
  entries := [];
  stats := rec(affine := 0, generic := 0, renamed := 0, perms := 0,
               ffe := 0, int := 0, failed := []);
  for nr in [1..n] do
    r := PRIMGRP_ReadLegacy(deg, nr);
    e := PRIMGRP_LegacyEntry(r, deg, nr);
    gen := PRIMGRP_GenericEntry(e, deg, nr);
    if gen <> fail then
      diff := Filtered([1..8], j -> e[j] <> gen[j]);
      if diff = [7] then
        # Field 7 only: the external data names A_n/S_n "A(n)" or leaves it
        # empty where the bundled data says "Alt(n)".  The reader overrides it
        # either way (lib/primitiv.gi:110), so this is dead data; normalise it.
        stats.renamed := stats.renamed + 1;
      elif diff <> [] then
        # Log field numbers, never values: field 2 of an A_n entry is a
        # 20000-digit integer.
        Add(stats.failed, [nr, "generic mismatch", diff]);
      fi;
      e := gen;
      stats.generic := stats.generic + 1;
    elif e[4] = "1" and Length(e[9]) > 0 and IsPerm(e[9][1]) then
      G := Group(e[9]);
      enum := PRIMGRP_AffineEnumeration(G, deg);
      rec_ := fail;
      if enum <> fail then
        rec_ := PRIMGRP_AffineMatrices(G, deg, enum);
      fi;
      if rec_ = fail or not PRIMGRP_CheckAffineRoundTrip(rec_, deg, enum, e[2]) then
        Add(stats.failed, [nr, "affine recovery"]);
        stats.perms := stats.perms + 1;
      else
        mats := rec_.mats;
        # "ffe" is the bare-matrix encoding data/gps*.g already uses; "int"
        # needs the tag, or the reader would rebuild a conjugate group.
        if enum = "ffe" then
          e[9] := mats;
        else
          e[9] := [ enum, mats ];
        fi;
        stats.affine := stats.affine + 1;
        stats.(enum) := stats.(enum) + 1;
      fi;
    else
      stats.perms := stats.perms + 1;
    fi;
    entries[nr] := e;
  od;
  return rec(entries := entries, stats := stats);
end);
