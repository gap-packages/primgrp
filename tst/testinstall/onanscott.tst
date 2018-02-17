gap> START_TEST("onanscott.tst");
gap> for i in [1..33] do
>      g:=PrimitiveGroup(3600,i);
>      h:=Group(GeneratorsOfGroup(g));
>      tg:=ONanScottType(g); th:=ONanScottType(h);
>      Print(i, " : ", tg, " : ", th, " : " , tg = th, "\n" );
>    od;
1 : 3b : 3b : true
2 : 3b : 3b : true
3 : 3b : 3b : true
4 : 3b : 3b : true
5 : 3b : 3b : true
6 : 4a : 4a : true
7 : 4b : 4b : true
8 : 4a : 4a : true
9 : 4a : 4a : true
10 : 4b : 4b : true
11 : 4b : 4b : true
12 : 4b : 4b : true
13 : 4b : 4b : true
14 : 4b : 4b : true
15 : 4b : 4b : true
16 : 4a : 4a : true
17 : 4b : 4b : true
18 : 4b : 4b : true
19 : 4b : 4b : true
20 : 4b : 4b : true
21 : 4b : 4b : true
22 : 4b : 4b : true
23 : 4b : 4b : true
24 : 4b : 4b : true
25 : 4b : 4b : true
26 : 4b : 4b : true
27 : 4b : 4b : true
28 : 4b : 4b : true
29 : 4b : 4b : true
30 : 4c : 4c : true
31 : 4c : 4c : true
32 : 4c : 4c : true
33 : 4c : 4c : true
gap> STOP_TEST( "onanscott.tst", 1);
