#############################################################################
##
##  dev/serialize.g
##
##  Compact writer for data/gps*.g.  Not read by read.g.
##
##  Two departures from GAP's pretty printer, both worth real bytes:
##  matrices over the prime field GF(p) are written as integer matrices scaled
##  once by Z(p)^0, and there is no indentation.  One entry per line, so data
##  fixes stay reviewable in a diff.
##

# A string literal, quoted and escaped.  Character 3 has to be spelled \c
# because that is how it got into the data in the first place: the name of
# PrimitiveGroup(625,657) is written "... 4 \circ 2^(1+4) ..." in gps4.g, and
# GAP's parser reads \c as character 3, leaving "4 <chr 3>irc".  Refuse to emit
# anything else unprintable rather than corrupt it silently.
BindGlobal("PRIMGRP_QuoteString", function(s)
  local out, c, i;
  out := "\"";
  for c in s do
    i := INT_CHAR(c);
    if c = '"' or c = '\\' then
      Add(out, '\\'); Add(out, c);
    elif c = '\n' then
      Append(out, "\\n");
    elif i = 3 then
      Append(out, "\\c");
    elif i < 32 or i > 126 then
      Error("cannot serialise character ", i, " in string \"", s, "\"");
    else
      Add(out, c);
    fi;
  od;
  Add(out, '"');
  return out;
end);

# String(obj) without GAP's spaces after commas and inside brackets.  Strings
# nested inside lists must come back out quoted -- the socle series in field 8
# is one, and an unquoted "Z" would read back as the Conway generator.
DeclareGlobalFunction("PRIMGRP_Compact");
InstallGlobalFunction(PRIMGRP_Compact, function(obj)
  if IsString(obj) and not IsEmpty(obj) then
    return PRIMGRP_QuoteString(obj);
  elif IsList(obj) then
    return Concatenation("[", JoinStringsWithSeparator(List(obj, PRIMGRP_Compact), ","), "]");
  fi;
  return String(obj);
end);

# A list of matrices over GF(p), as  Z(p)^0*[[..],[..]]  with integer entries.
BindGlobal("PRIMGRP_CompactMatrices", function(mats, p)
  if Length(mats) = 0 then
    return "[]";
  fi;
  return Concatenation("Z(", String(p), ")^0*",
                       PRIMGRP_Compact(List(mats, m -> List(m, r -> List(r, IntFFE)))));
end);

BindGlobal("PRIMGRP_CompactAffine", function(mats, deg)
  local fac;
  fac := Factors(deg);
  if Length(mats) > 0 and IsFFE(mats[1][1][1]) then
    return PRIMGRP_CompactMatrices(mats, fac[1]);
  fi;
  # already integer matrices, e.g. straight from dev/affine.g
  return Concatenation("Z(", String(fac[1]), ")^0*", PRIMGRP_Compact(mats));
end);

# Field 9 of an entry: markers, matrices, permutations, or nothing.
BindGlobal("PRIMGRP_CompactGenerators", function(gens, deg, onanscott)
  if IsString(gens) and not IsEmpty(gens) then  # "Alt", "Sym", "psl", "pgl"
    return PRIMGRP_QuoteString(gens);
  elif Length(gens) = 0 then
    return "[]";
  elif IsString(gens[1]) then
    # ["int", mats]: affine, but enumerating F_p^d by base-p digit value rather
    # than by GAP's order on field elements.  See dev/affine.g.
    return Concatenation("[\"", gens[1], "\",",
                         PRIMGRP_CompactAffine(gens[2], deg), "]");
  elif not IsPerm(gens[1]) then                # affine: matrices over GF(p)
    return PRIMGRP_CompactAffine(gens, deg);
  fi;
  return PRIMGRP_Compact(gens);                # permutations
end);

# One PRIMGRP[deg][nr] entry.  <entry> is the 9- or 10-element list.
#
# Groups we know generically collapse to PGAlt(deg,nr) and friends -- but only
# when the stored entry agrees with the generic rule in every field.  Five
# entries do not (degrees 2, 3, 4, where A_n/S_n really is affine and the socle
# really is elementary abelian) and they are written out in full.  Silently
# generalising over a disagreement is how wrong data gets baked in.
BindGlobal("PRIMGRP_CompactEntry", function(entry, deg)
  local parts, generic, tag;
  if IsBound(PRIMGRP_GenericEntry) and not IsBound(entry[10]) then
    generic := PRIMGRP_GenericEntry(entry, deg, entry[1]);
    if generic <> fail and generic = entry then
      tag := rec(Alt := "PGAlt", Sym := "PGSym",
                 psl := "PGPsl", pgl := "PGPgl").(entry[9]);
      return Concatenation(tag, "(", String(deg), ",", String(entry[1]), ")");
    fi;
  fi;
  parts := [ String(entry[1]),
             String(entry[2]),
             String(entry[3]),
             PRIMGRP_QuoteString(entry[4]),
             PRIMGRP_Compact(entry[5]),
             String(entry[6]),
             PRIMGRP_QuoteString(entry[7]),
             PRIMGRP_Compact(entry[8]),
             PRIMGRP_CompactGenerators(entry[9], deg, entry[4]) ];
  if IsBound(entry[10]) then
    Add(parts, String(entry[10]));
  fi;
  return Concatenation("[", JoinStringsWithSeparator(parts, ","), "]");
end);

# The whole PRIMGRP[deg]:=[...]; assignment, one entry per line.
BindGlobal("PRIMGRP_CompactDegree", function(entries, deg)
  local out, i;
  out := [ "PRIMGRP[", String(deg), "]:=[\n" ];
  for i in [1..Length(entries)] do
    Add(out, PRIMGRP_CompactEntry(entries[i], deg));
    if i < Length(entries) then Add(out, ",\n"); fi;
  od;
  Add(out, "];\n");
  return Concatenation(out);
end);
