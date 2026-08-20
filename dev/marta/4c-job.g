#############################################################################
##
##  dev/marta/4c-job.g
##
##  One MARTA job = one type 4c entry.  Needs marta_deg and marta_nr via -c.
##
##  The result is the wreath decomposition as a flat integer list, which is
##  what MARTA summarises natively:
##
##      m, k, ngens, then per generator the k images lists of length m
##      followed by sigma as a list of length k.
##
##  An entry that cannot be decomposed writes [ 0 ] and exits 0 -- a refusal,
##  not a failure, so MARTA does not retry it.
##
EXITSTATUS := 0;;
if not IsBound(marta_deg) or not IsBound(marta_nr) then
    Print("marta_deg and marta_nr must be set via GAP -c\n");  QuitGap(1);
fi;

Read("dev/product4c.g");;
BreakOnError := false;;
start_time := Runtime();;

res := CALL_WITH_CATCH(PRIMGRP_Extract4c,
         [PrimitiveGroup(marta_deg, marta_nr), marta_deg,
          PRIMGrp(marta_deg, marta_nr)[8][3]]);;
flat := [0];;
if res[1] and not IsBound(res[2].err) then
    r := res[2];;
    # check the rebuild is the original relabelled before reporting anything
    G := PrimitiveGroup(marta_deg, marta_nr);;
    if PGProductAction4c(r.m, r.k, r.els)
       = Image(ConjugatorIsomorphism(G, r.rel), G) then
        flat := [r.m, r.k, Length(r.els)];;
        for e in r.els do
            for i in [1..r.k] do
                Append(flat, OnTuples([1..r.m], e[i]));
            od;
            Append(flat, OnTuples([1..r.k], e[r.k+1]));
        od;
    else
        Print("REFUSED rebuild differs from the relabelled original\n");
    fi;
elif res[1] then
    Print("REFUSED ", res[2].err, "\n");
else
    Print("REFUSED exception during extraction\n");
fi;

output := OutputTextFile(GAPInfo.SystemEnvironment.MARTA_RESULT_PATH, false);;
for value in flat do AppendTo(output, String(value), "\n"); od;
CloseStream(output);

metadata := OutputTextFile(GAPInfo.SystemEnvironment.MARTA_METADATA_PATH, false);;
AppendTo(metadata, "{\n");
AppendTo(metadata, "  \"software_versions\": [\n");
AppendTo(metadata, "    {\"name\": \"GAP\", \"version\": \"", GAPInfo.Version, "\", \"metadata\": {}}\n");
AppendTo(metadata, "  ],\n");
AppendTo(metadata, "  \"timings\": {\"inner_runtime_seconds\": ", String(Runtime()-start_time), "e-3}\n");
AppendTo(metadata, "}\n");
CloseStream(metadata);
QuitGap(EXITSTATUS);
