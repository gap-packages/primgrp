#############################################################################
##
##  dev/generic.g
##
##  Entries for groups we know generically, so the data files stop repeating
##  what is determined by the degree alone (issue #53).  These build the full
##  9-element PRIMGRP entry at read time, which keeps the fields that
##  PrimitiveGroupsIterator and PrimitiveIdentification read directly out of
##  PRIMGRP available without constructing the group.
##
##  Destined for lib/primitiv.gi; the short names are deliberate, they occur
##  ~9300 times across the data files.
##

BindGlobal("PGAlt", function(deg, nr)
  return [ nr, Factorial(deg)/2, 1, "2", [[deg-1,1]], deg-2,
           Concatenation("Alt(", String(deg), ")"), ["A",deg,1], "Alt" ];
end);

BindGlobal("PGSym", function(deg, nr)
  return [ nr, Factorial(deg), 0, "2", [[deg-1,1]], deg,
           Concatenation("Sym(", String(deg), ")"), ["A",deg,1], "Sym" ];
end);

# Degree q+1, natural action on the projective line.
BindGlobal("PGPsl", function(deg, nr)
  local q;
  q := deg-1;
  return [ nr, q*(q^2-1)/Gcd(2,q-1), 1, "2", [[q,1]], 2,
           Concatenation("PSL(2, ", String(q), ")"), ["L",[2,q],1], "psl" ];
end);

BindGlobal("PGPgl", function(deg, nr)
  local q;
  q := deg-1;
  return [ nr, q*(q^2-1), 0, "2", [[q,1]], 3,
           Concatenation("PGL(2, ", String(q), ")"), ["L",[2,q],1], "pgl" ];
end);

BindGlobal("PRIMGRP_GenericEntry", function(entry, deg, nr)
  if entry[9] = "Alt" then return PGAlt(deg, nr);
  elif entry[9] = "Sym" then return PGSym(deg, nr);
  elif entry[9] = "psl" then return PGPsl(deg, nr);
  elif entry[9] = "pgl" then return PGPgl(deg, nr);
  fi;
  return fail;
end);
