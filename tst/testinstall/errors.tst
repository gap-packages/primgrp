#@local it
gap> START_TEST("errors.tst");

# Invalid arguments are refused with a message naming the argument, not with
# an error from deep inside a list access.

#
gap> PrimitiveGroup(0, 1);
Error, PrimitiveGroup: <deg> must be a positive integer
gap> PrimitiveGroup(-3, 1);
Error, PrimitiveGroup: <deg> must be a positive integer
gap> PrimitiveGroup("5", 1);
Error, PrimitiveGroup: <deg> must be a positive integer
gap> PrimitiveGroup(8192, 1);
Error, PrimitiveGroup: <deg> must be at most 8191
gap> PrimitiveGroup(5, 0);
Error, PrimitiveGroup: <nr> must be a positive integer
gap> PrimitiveGroup(5, 1/2);
Error, PrimitiveGroup: <nr> must be a positive integer
gap> PrimitiveGroup(5, 6);
Error, PrimitiveGroup: there are only 5 primitive groups of degree 5
gap> PrimitiveGroup(1, 1);
Error, PrimitiveGroup: there are only 0 primitive groups of degree 1

#
gap> NrPrimitiveGroups(0);
Error, NrPrimitiveGroups: <deg> must be a positive integer
gap> NrPrimitiveGroups(-1);
Error, NrPrimitiveGroups: <deg> must be a positive integer
gap> NrPrimitiveGroups("x");
Error, NrPrimitiveGroups: <deg> must be a positive integer
gap> NrPrimitiveGroups(8192);
Error, NrPrimitiveGroups: <deg> must be at most 8191

#
gap> PrimitiveIdentification(Group(()));
Error, PrimitiveIdentification: <G> must move at least two points
gap> PrimitiveIdentification(Group((1,2), (3,4)));
Error, Group must operate primitively
gap> PrimitiveIdentification(Group((1,2)(3,4), (1,3)(2,4)));
Error, Group must operate primitively
gap> SimsNo(Group((1,2)(3,4), (1,3)(2,4)));
Error, Group must operate primitively
gap> SimsNo(SymmetricGroup(51));
Error, SimsNo is defined only for primitive groups of degree <= 50

# The selection functions share the argument parsing of the iterator.
gap> AllPrimitiveGroups(NrMovedPoints);
Error, the arguments must be pairs of a function and a value
gap> it:= PrimitiveGroupsIterator(NrMovedPoints, 5, Size);
Error, the arguments must be pairs of a function and a value
gap> AllPrimitiveGroups(5, IsSolvable);
Error, argument 1 must be a function
gap> OnePrimitiveGroup(NrMovedPoints, 5, 60, Size);
Error, argument 3 must be a function
gap> AllPrimitiveGroups(NrMovedPoints, "x");
Error, the value for NrMovedPoints must be integers or a function
gap> AllPrimitiveGroups(NrMovedPoints, 5, Size, "60");
Error, the value for Size must be integers or a function

#
gap> STOP_TEST("errors.tst");
