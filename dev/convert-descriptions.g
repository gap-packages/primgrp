#############################################################################
##
##  dev/convert-descriptions.g
##
##  Name the construction of an entry instead of calling it.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-descriptions.g
##
##  An entry built by a construction was written as the construction itself,
##  PGAlt, which GAP evaluates while reading the file.  What lands in PRIMGRP
##  is then a function, and a pass that writes the entries back out has only
##  NameFunction to describe it.  That works for a bare name and fails for a
##  call: PGAltOnSets(10,2) comes back as "unknown", its arguments gone.
##
##  So the file says ["Alt"], and ["PSL",2,5] where it said PGPsl(2,5).  The
##  arguments now have somewhere to live, which is the point; PRIMGrp matches
##  the name against the constructions the library offers, so a data file can
##  ask for one of those and for nothing else.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

##  The constructions a data file may name, and what it used to write instead.
PRIMGRP_Bare := rec(PGAlt := "Alt", PGSym := "Sym");

PRIMGRP_Called := rec(PGAltOnSets := "AltOnSets", PGSymOnSets := "SymOnSets",
                      PGPsl := "PSL", PGPgl := "PGL",
                      PGPsigmaL := "PSigmaL", PGPgammaL := "PGammaL");

##  Turn one line into a description, or return fail.
PRIMGRP_Describe := function(line)
  local body, open, name, args;
  if Length(line) < 2 or line[Length(line)] <> ',' then
    return fail;
  fi;
  body := line{[1..Length(line)-1]};

  if body in RecNames(PRIMGRP_Bare) then
    return Concatenation("[\"", PRIMGRP_Bare.(body), "\"]");
  fi;

  open := Position(body, '\(');
  if open = fail or body[Length(body)] <> ')' then
    return fail;
  fi;
  name := body{[1..open-1]};
  if not name in RecNames(PRIMGRP_Called) then
    return fail;
  fi;
  args := body{[open+1..Length(body)-1]};
  return Concatenation("[\"", PRIMGRP_Called.(name), "\",", args, "]");
end;

##  A content line is one ending in a comma that begins with neither
##  "PRIMGRP[" nor "];".  That is what the normalisation pass bought.
PRIMGRP_ConvertFile := function(path)
  local lines, out, line, text, n;
  lines := SplitString(StringFile(path), "\n");
  out := [];
  n := 0;
  for line in lines do
    if Length(line) = 0 then
      continue;
    fi;
    text := PRIMGRP_Describe(line);
    if text = fail then
      Add(out, line);
    else
      Add(out, Concatenation(text, ","));
      n := n + 1;
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
