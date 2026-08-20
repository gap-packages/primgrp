#############################################################################
##
##  dev/marta/sweep-job.g
##
##  One MARTA job = one degree of the consistency sweep.  Needs marta_deg set
##  via -c, and the primgrp worktree reachable as a GAP root (-l).
##
##  Result is [ degree, groups, complaints ]; the complaints themselves go to
##  stdout, which MARTA keeps with the attempt.
##
EXITSTATUS := 0;;
if not IsBound(marta_deg) then
    Print("marta_deg must be set via GAP -c\n");  QuitGap(1);
fi;
if not IsBound(GAPInfo.SystemEnvironment.MARTA_RESULT_PATH) then
    Print("MARTA_RESULT_PATH must be set by the environment\n");  QuitGap(1);
fi;

LoadPackage("primgrp");;
ReadPackage("primgrp", "tst/testutils.g");;
SetInfoLevel(InfoWarning, 0);;

start_time := Runtime();;
bad := PrimGrpCheckDegree(marta_deg, PrimGrpCheckAll);;
for x in bad do Print("COMPLAINT ", x, "\n"); od;

output := OutputTextFile(GAPInfo.SystemEnvironment.MARTA_RESULT_PATH, false);;
for value in [marta_deg, NrPrimitiveGroups(marta_deg), Length(bad)] do
    AppendTo(output, String(value), "\n");
od;
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
