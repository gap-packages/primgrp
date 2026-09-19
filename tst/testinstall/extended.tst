gap> START_TEST("extended.tst");

# Degrees 4096 to 8191 ship with the package.
gap> ForAll([4096, 8191], PrimitiveGroupsAvailable);
true
gap> Sum([4096..8191], NrPrimitiveGroups);
29776

# Affine by matrices and diagonal type by generators, their orders computed
# from the generators rather than set from the entry.  Not so for prime degree:
# the translation is one 8191-cycle, and the stabiliser chain takes minutes.
# 4913 = 17^3 is read in the archive's numbering of GF(17)^3 and built in
# GAP's, so the group is a conjugate of the archive's, of the same order.
gap> List([[4096,1],[4913,1],[5616,1]],
>         p -> Size(Group(GeneratorsOfGroup(PrimitiveGroup(p[1],p[2])))));
[ 53248, 58956, 31539456 ]
gap> PRIMGrp(8191,48){[2,7]};
[ 67084290, "AGL(1, 8191)" ]

# Entries written as descriptions: a product action from its wreath elements,
# the L series on points, PSL(2,4096) extended by the 6th power of the
# Frobenius, Alt(92) on pairs, PGL(2,97) on pairs, and PSL(4,8) on lines.
gap> PRIMGrp(4096,1392)[9][1];
"4c"
gap> PRIMGrp(4097,5){[7,9]};
[ "PSL(2,4096)", "psl" ]
gap> PRIMGrp(4097,6)[9];
[ "PSL", 2, 4096, [ [ 0, 6 ] ], "" ]
gap> PRIMGrp(4186,1){[7,9]};
[ "A(92)", [ "sets", [ "Alt", 92 ], 2 ] ]
gap> PRIMGrp(4745,1){[7,9]};
[ "PSL(4,8)", [ "subspaces", [ "PSL", 4, 8 ], 2 ] ]
gap> ForAll([[4096,1392],[4097,6],[4753,2],[4745,1]], p ->
>      Size(Group(GeneratorsOfGroup(PrimitiveGroup(p[1],p[2]))))
>      = PRIMGrp(p[1],p[2])[2]);
true

#
gap> STOP_TEST("extended.tst", 1);
