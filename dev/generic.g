#############################################################################
##
##  dev/generic.g
##
##  What the writer needs: given a stored entry, is it one the degree (and, for
##  the projective groups, the dimension) determines?  It builds the expected
##  entry itself rather than calling PGAlt and friends from lib/primitiv.gi,
##  so that regeneration works whichever form the data currently on disk uses.
##

PRIMGRP_PslOrder := function(dim, q)
  local o, i;
  o := q^(dim*(dim-1)/2);
  for i in [2..dim] do o := o * (q^i-1); od;
  return o / Gcd(dim, q-1);
end;

##  Fields 1..8 of the entry the generic form would produce, or fail.
PRIMGRP_GenericFields := function(entry, deg, nr)
  local q, t, d;
  if entry[9] = "Alt" then
    return [ nr, Factorial(deg)/2, 1, "2", [[deg-1,1]], deg-2,
             Concatenation("Alt(", String(deg), ")"), ["A",deg,1] ];
  elif entry[9] = "Sym" then
    return [ nr, Factorial(deg), 0, "2", [[deg-1,1]], deg,
             Concatenation("Sym(", String(deg), ")"), ["A",deg,1] ];
  elif entry[9] = "psl" or (IsList(entry[9]) and Length(entry[9]) = 3
                            and entry[9][1] = "psl") then
    q := deg-1;
    if not IsString(entry[9]) then q := entry[9][3]; fi;
    t := 2;
    if q mod 2 = 0 then t := 3; fi;
    d := 2;
    if not IsString(entry[9]) then d := entry[9][2]; fi;
    if d = 2 and q mod 2 = 0 then t := 3; else t := 2; fi;
    return [ nr, PRIMGRP_PslOrder(d,q), 1, "2", [[deg-1,1]], t,
             Concatenation("PSL(", String(d), ", ", String(q), ")"), ["L",[d,q],1] ];
  elif entry[9] = "pgl" or (IsList(entry[9]) and Length(entry[9]) = 3
                            and entry[9][1] = "pgl") then
    q := deg-1;
    if not IsString(entry[9]) then q := entry[9][3]; fi;
    d := 2;
    if not IsString(entry[9]) then d := entry[9][2]; fi;
    return [ nr, PRIMGRP_PslOrder(d,q)*Gcd(d,q-1), 0, "2", [[deg-1,1]], 3,
             Concatenation("PGL(", String(d), ", ", String(q), ")"), ["L",[d,q],1] ];
  fi;
  return fail;
end;
