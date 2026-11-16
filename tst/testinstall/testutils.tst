gap> START_TEST("primutils.tst");

#
# Disable warnings which depend on Conway Polynomial databases 
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
gap> ReadPackage("primgrp", "lib/testutils.g");
true
gap> for d in [1,2,3,4] do
>   for q in [2,3,4,5,7,8,9] do
>     facs := Factors(q);
>     p := facs[1];
>     e := Length(facs);
>     Assert(0, Size(SL(d,q)) * q^d = Size(ASL(d,q)));
>     Assert(0, Size(GL(d,q)) * q^d = Size(AGL(d,q)));
>     Assert(0, Size(SigmaL(d,q)) * q^d = Size(ASigmaL(d,q)));
>     Assert(0, Size(GammaL(d,q)) * q^d = Size(AGammaL(d,q)));
>     if d > 1 then
>         Assert(0, Size(PSL(d,q)) * e = Size(PSigmaL(d,q)));
>         Assert(0, Size(PGL(d,q)) * e = Size(PGammaL(d,q)));
>     fi;
>   od;
> od;
gap> 

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "primutils.tst", 1);
