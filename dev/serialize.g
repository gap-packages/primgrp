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

##  One d x d matrix over GF(p), read row by row as a base-p number, written
##  in base 256 and base64 encoded.  Cuts an 8 x 8 matrix over GF(3) from
##  about 150 characters to 22.
BindGlobal("PRIMGRP_PackMatrix", function(p, d, m)
  local n, r, x, len, bytes, k, i;
  n := 0;
  for r in m do
    for x in r do
      n := n*p + x;
    od;
  od;
  len := 1;
  while p^(d*d) > 256^len do len := len+1; od;
  bytes := "";
  k := n;
  for i in [1..len] do
    bytes := Concatenation([CHAR_INT(k mod 256)], bytes);
    k := QuoInt(k, 256);
  od;
  return Base64String(bytes);
end);

BindGlobal("PRIMGRP_CompactAffine", function(mats, deg)
  local fac, p, d, ints;
  fac := Factors(deg);
  p := fac[1]; d := Length(fac);
  if Length(mats) = 0 then
    return "[]";
  fi;
  if IsFFE(mats[1][1][1]) then
    ints := List(mats, m -> List(m, r -> List(r, IntFFE)));
  else
    ints := mats;
  fi;
  # packed unless the matrices are so small that the packing is not a saving
  if d*d >= 9 then
    return Concatenation("[\"b64\",\"ffe\",",
             PRIMGRP_Compact(List(ints, m -> PRIMGRP_PackMatrix(p, d, m))), "]");
  fi;
  return Concatenation("Z(", String(p), ")^0*", PRIMGRP_Compact(ints));
end);

# Field 9 of an entry: markers, matrices, permutations, or nothing.
BindGlobal("PRIMGRP_CompactGenerators", function(gens, deg, onanscott)
  local fac;
  if IsString(gens) and not IsEmpty(gens) then  # "Alt", "Sym", "psl", "pgl"
    return PRIMGRP_QuoteString(gens);
  elif Length(gens) = 0 then
    return "[]";
  elif IsList(gens) and Length(gens) = 3 and IsString(gens[1])
       and (gens[1] = "psl" or gens[1] = "pgl") then
    return PRIMGRP_Compact(gens);          # ["psl",dim,q], written as data
  elif IsString(gens[1]) then
    # ["int", mats]: affine, but enumerating F_p^d by base-p digit value rather
    # than by GAP's order on field elements.  See dev/affine.g.
    # ["int", mats]: pack the same way, recording the enumeration
    fac := Factors(deg);
    if Length(fac)^2 >= 9 and Length(gens[2]) > 0 then
      return Concatenation("[\"b64\",\"", gens[1], "\",",
               PRIMGRP_Compact(List(gens[2],
                 m -> PRIMGRP_PackMatrix(fac[1], Length(fac),
                        List(m, r -> List(r, function(x)
                          if IsFFE(x) then return IntFFE(x); fi; return x; end))))), "]");
    fi;
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
  local parts, generic;
  # PRIMGRP_GenericFields and PRIMGRP_PglIdentified live in dev/generic.g and
  # are optional here.  Reached through ValueGlobal so that reading this file
  # without them is silent rather than a parse-time warning.
  if IsBoundGlobal("PRIMGRP_GenericFields") then
    generic := ValueGlobal("PRIMGRP_GenericFields")(entry, deg, entry[1]);
    # Fields 1..8 must agree exactly.  Field 9 need not: it is what the
    # constructor replaces, whether that was the marker "psl" or a list of
    # permutations spelling out PSL(d,q) the long way.
    # Fields 1..6 and 8 must agree exactly.  Field 7, the name, may be filled
    # in where the entry has none -- 497 of these groups are unnamed, and a
    # correct name is better than no name -- but never overwritten.
    if generic <> fail
       and generic.fields{[1..6]} = entry{[1..6]}
       and generic.fields[8] = entry[8]
       and (generic.fields[7] = entry[7] or entry[7] = ""
            or (IsBoundGlobal("PRIMGRP_PglIdentified")
                and [deg,entry[1]] in ValueGlobal("PRIMGRP_PglIdentified"))) then
      return generic.text;
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
  # Field 10, Sims' number, has moved to PRIMGRP_SIMSNO in lib/primitiv.grp.
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
