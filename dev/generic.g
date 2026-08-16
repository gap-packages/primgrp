#############################################################################
##
##  dev/generic.g
##
##  PGAlt, PGSym, PGPsl and PGPgl live in lib/primitiv.gi -- the data files
##  call them, so they have to be there.  What belongs here is only the
##  question the writer asks: is this entry one of them?
##

BindGlobal("PRIMGRP_GenericEntry", function(entry, deg, nr)
  if entry[9] = "Alt" then return PGAlt(deg, nr);
  elif entry[9] = "Sym" then return PGSym(deg, nr);
  elif entry[9] = "psl" then return PGPsl(deg, nr);
  elif entry[9] = "pgl" then return PGPgl(deg, nr);
  fi;
  return fail;
end);
