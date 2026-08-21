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
  local q, t, d, e, sg;
  if entry[9] = "Alt" then
    return rec(text := "PGAlt", fields :=
           [ nr, Factorial(deg)/2, 1, "2", [[deg-1,1]], deg-2,
             Concatenation("Alt(", String(deg), ")"), ["A",deg,1] ]);
  elif entry[9] = "Sym" then
    return rec(text := "PGSym", fields :=
           [ nr, Factorial(deg), 0, "2", [[deg-1,1]], deg,
             Concatenation("Sym(", String(deg), ")"), ["A",deg,1] ]);
  elif entry[9] = "psl" or (IsList(entry[9]) and Length(entry[9]) = 3
                            and entry[9][1] = "psl") then
    q := deg-1;
    if not IsString(entry[9]) then q := entry[9][3]; fi;
    t := 2;
    if q mod 2 = 0 then t := 3; fi;
    d := 2;
    if not IsString(entry[9]) then d := entry[9][2]; fi;
    if d = 2 and q mod 2 = 0 then t := 3; else t := 2; fi;
    return rec(text := Concatenation("PGPsl(", String(d), ",", String(q), ")"),
               fields := [ nr, PRIMGRP_PslOrder(d,q), 1, "2", [[deg-1,1]], t,
             Concatenation("PSL(", String(d), ", ", String(q), ")"), ["L",[d,q],1] ]);
  elif entry[9] = "pgl" or (IsList(entry[9]) and Length(entry[9]) = 3
                            and entry[9][1] = "pgl") then
    q := deg-1;
    if not IsString(entry[9]) then q := entry[9][3]; fi;
    d := 2;
    if not IsString(entry[9]) then d := entry[9][2]; fi;
    return rec(text := Concatenation("PGPgl(", String(d), ",", String(q), ")"),
               fields := [ nr, PRIMGRP_PslOrder(d,q)*Gcd(d,q-1), 0, "2", [[deg-1,1]], 3,
             Concatenation("PGL(", String(d), ", ", String(q), ")"), ["L",[d,q],1] ]);
  fi;
  # PSigmaL and PGammaL, likewise only on an explicit list: their orders are
  # |PSL|*e and |PGL|*e for q = p^e, and at many degrees a sibling extension
  # has the same order.
  sg := fail;
  if IsBoundGlobal("PRIMGRP_SigmaGamma") then
    sg := ValueGlobal("PRIMGRP_SigmaGamma");
  fi;
  if sg <> fail and IsBound(sg.(String(deg)))
     and IsBound(sg.(String(deg)).(String(nr)))
     and IsList(entry[8]) and entry[8][1] = "L" and IsList(entry[8][2])
     and entry[8][3] = 1 and not IsString(entry[9]) then
    d := entry[8][2][1]; q := entry[8][2][2];
    e := Length(Factors(q));
    if sg.(String(deg)).(String(nr)) = "sigma" then
      return rec(text := Concatenation("PGPsigmaL(", String(d), ",", String(q), ")"),
                 fields := [ nr, PRIMGRP_PslOrder(d,q)*e, 0, "2", [[deg-1,1]], 2,
                             Concatenation("PSigmaL(", String(d), ", ", String(q), ")"),
                             ["L",[d,q],1] ]);
    else
      t := 2;
      if d = 2 then t := 3; fi;
      return rec(text := Concatenation("PGPgammaL(", String(d), ",", String(q), ")"),
                 fields := [ nr, PRIMGRP_PslOrder(d,q)*Gcd(d,q-1)*e, 0, "2",
                             [[deg-1,1]], t,
                             Concatenation("PGammaL(", String(d), ", ", String(q), ")"),
                             ["L",[d,q],1] ]);
    fi;
  fi;

  # A permutation-stored group whose socle is L(d,q), on the (q^d-1)/(q-1)
  # points of PG(d-1,q), and whose order is |PGL(d,q)|, is PGL(d,q) -- but
  # only where no other entry of that degree has the same socle and order.
  # At degree 170, PGL(2,169), PSigmaL(2,169) and PSL(2,169).2_3 all do, and
  # they are different groups; PRIMGRP_PglUnambiguous lists the safe ones.
  # The three tables below are optional side inputs; reaching them through
  # ValueGlobal keeps reading this file silent when they are absent.
  if ((IsBoundGlobal("PRIMGRP_PglUnambiguous")
       and [deg,nr] in ValueGlobal("PRIMGRP_PglUnambiguous"))
      or (IsBoundGlobal("PRIMGRP_PglIdentified")
          and [deg,nr] in ValueGlobal("PRIMGRP_PglIdentified")))
     and IsList(entry[8]) and entry[8][1] = "L" and IsList(entry[8][2])
     and entry[8][3] = 1 and not IsString(entry[9]) then
    d := entry[8][2][1]; q := entry[8][2][2];
    if deg = (q^d-1)/(q-1)
       and entry[2] = PRIMGRP_PslOrder(d,q)*Gcd(d,q-1) then
      return rec(text := Concatenation("PGPgl(", String(d), ",", String(q), ")"),
                 fields := [ nr, PRIMGRP_PslOrder(d,q)*Gcd(d,q-1), 0, "2",
                             [[deg-1,1]], 3,
                             Concatenation("PGL(", String(d), ", ", String(q), ")"),
                             ["L",[d,q],1] ]);
    fi;
  fi;

  # A permutation-stored group whose socle is L(d,q), on the (q^d-1)/(q-1)
  # points of PG(d-1,q), and whose order is |PSL(d,q)|, is PSL(d,q) itself.
  # It is conjugate to GAP's PSL(d,q) because no other entry of that degree
  # shares this socle and order, and the library is complete.
  if IsList(entry[8]) and entry[8][1] = "L" and IsList(entry[8][2])
     and entry[8][3] = 1 and not IsString(entry[9]) then
    d := entry[8][2][1]; q := entry[8][2][2];
    if deg = (q^d-1)/(q-1) and entry[2] = PRIMGRP_PslOrder(d,q) then
      t := 2;
      if d = 2 and q mod 2 = 0 then t := 3; fi;
      return rec(text := Concatenation("PGPsl(", String(d), ",", String(q), ")"),
                 fields := [ nr, PRIMGRP_PslOrder(d,q), 1, "2", [[deg-1,1]], t,
               Concatenation("PSL(", String(d), ", ", String(q), ")"),
               ["L",[d,q],1] ]);
    fi;
  fi;
  return fail;
end;
