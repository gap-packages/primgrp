#############################################################################
##
##  dev/probe-compress.g
##
##  How much is there to gain from the issue #56 encoding, and at what cost?
##
##  Only ~4600 of the 54334 groups are still stored as explicit permutations,
##  but they are essentially all of the remaining bulk.  This tries the small
##  degree representation on them, biggest entry first, so the question "is
##  phase 4 worth it" is answered by the early results rather than the last.
##
##  Resumable, and with an effective time budget: a "start" marker is written
##  before each attempt and a "done" marker after.  A group whose process was
##  killed mid-attempt therefore has a start and no done, and is skipped on the
##  next run rather than hanging it again.  The budget is the outer timeout.
##

LoadPackage("primgrp");

##  Every group still stored as permutations, heaviest first.
BindGlobal("PRIMGRP_PermStoredGroups", function(degrees)
  local out, deg, nr, l;
  out := [];
  for deg in degrees do
    PrimGrpLoad(deg);
    for nr in [1..NrPrimitiveGroups(deg)] do
      l := PRIMGrp(deg, nr);
      if not IsString(l[9]) and Length(l[9]) > 0 and IsPerm(l[9][1]) then
        Add(out, [Length(String(l[9])), deg, nr]);
      fi;
    od;
  od;
  Sort(out, function(a, b) return a[1] > b[1]; end);
  return out;
end);

##  A word as one character per generator, upper case for an inverse, as in
##  issue #56: x^-4*y^-1*x -> "AAAABa".  GAP's own printing of a free group
##  element is several times longer, and the words are most of the payload.
BindGlobal("PRIMGRP_EncodeWord", function(w)
  local ext, res, i, c, e;
  ext := ExtRepOfObj(w);
  res := "";
  for i in [1, 3 .. Length(ext)-1] do
    c := CharInt(IntChar('a') + ext[i] - 1);
    e := ext[i+1];
    if e < 0 then
      c := UppercaseChar(c);
      e := -e;
    fi;
    Append(res, ListWithIdenticalEntries(e, c));
  od;
  return res;
end);

##  Iterated SmallerDegreePermutationRepresentation, as in issue #56, plus the
##  point stabiliser as words.  Returns the encoded size, or fail.
BindGlobal("PRIMGRP_TryCompress", function(G, rounds)
  local iso, i, H, S, hom, words, w;
  iso := IdentityMapping(G);
  for i in [1..rounds] do
    iso := iso * SmallerDegreePermutationRepresentation(Image(iso));
  od;
  H := Image(iso);
  if NrMovedPoints(H) >= NrMovedPoints(G) then
    return fail;
  fi;
  S := Image(iso, Stabilizer(G, 1));
  hom := EpimorphismFromFreeGroup(H);
  words := List(GeneratorsOfGroup(S), x -> PreImagesRepresentative(hom, x));
  return rec(degree := NrMovedPoints(H),
             bytes := Length(String(GeneratorsOfGroup(H)))
                      + Sum(words, w -> Length(PRIMGRP_EncodeWord(w)) + 3),
             rawbytes := Length(String(GeneratorsOfGroup(H))) + Length(String(words)));
end);

BindGlobal("PRIMGRP_Probe", function(items, markDir, logFile)
  local it, deg, nr, before, G, t, r, mark, out, line;
  for it in items do
    before := it[1]; deg := it[2]; nr := it[3];
    mark := Concatenation(markDir, "/g", String(deg), "_", String(nr));
    if IsReadableFile(Concatenation(mark, ".start")) then
      continue;                      # attempted before; done or too slow
    fi;
    FileString(Concatenation(mark, ".start"), "");
    G := PrimitiveGroup(deg, nr);
    t := Runtime();
    r := PRIMGRP_TryCompress(G, 5);
    if r = fail then
      line := Concatenation(String(deg), " ", String(nr), " before=", String(before),
                            " NOGAIN ms=", String(Runtime()-t), "\n");
    else
      line := Concatenation(String(deg), " ", String(nr), " before=", String(before),
                            " after=", String(r.bytes), " raw=", String(r.rawbytes),
                            " smalldeg=", String(r.degree),
                            " ms=", String(Runtime()-t), "\n");
    fi;
    out := OutputTextFile(logFile, true);
    SetPrintFormattingStatus(out, false);
    WriteAll(out, line);
    CloseStream(out);
    FileString(Concatenation(mark, ".done"), "");
  od;
end);

PRIMGRP_Probe(PRIMGRP_PermStoredGroups(PRIMGRP_Degrees),
              PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
