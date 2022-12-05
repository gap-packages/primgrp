gap> START_TEST("consistency99.tst");

#
# Disable warnings which depend on Conway Polynomial databases 
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
# Define a function to check the primitive groups of degree n
#
gap> ReadPackage("primgrp", "lib/testutils.g");
true
gap> PrimGrpConsistencyCheckDegree(2);
gap> PrimGrpConsistencyCheckDegree(3);
gap> PrimGrpConsistencyCheckDegree(4);
gap> PrimGrpConsistencyCheckDegree(5);
gap> PrimGrpConsistencyCheckDegree(6);
gap> PrimGrpConsistencyCheckDegree(7);
gap> PrimGrpConsistencyCheckDegree(8);
gap> PrimGrpConsistencyCheckDegree(9);
gap> PrimGrpConsistencyCheckDegree(10);
gap> PrimGrpConsistencyCheckDegree(11);
gap> PrimGrpConsistencyCheckDegree(12);
gap> PrimGrpConsistencyCheckDegree(13);
gap> PrimGrpConsistencyCheckDegree(14);
gap> PrimGrpConsistencyCheckDegree(15);
gap> PrimGrpConsistencyCheckDegree(16);
gap> PrimGrpConsistencyCheckDegree(17);
gap> PrimGrpConsistencyCheckDegree(18);
gap> PrimGrpConsistencyCheckDegree(19);
gap> PrimGrpConsistencyCheckDegree(20);
gap> PrimGrpConsistencyCheckDegree(21);
gap> PrimGrpConsistencyCheckDegree(22);
gap> PrimGrpConsistencyCheckDegree(23);
gap> PrimGrpConsistencyCheckDegree(24);
gap> PrimGrpConsistencyCheckDegree(25);
gap> PrimGrpConsistencyCheckDegree(26);
gap> PrimGrpConsistencyCheckDegree(27);
gap> PrimGrpConsistencyCheckDegree(28);
gap> PrimGrpConsistencyCheckDegree(29);
gap> PrimGrpConsistencyCheckDegree(30);
gap> PrimGrpConsistencyCheckDegree(31);
gap> PrimGrpConsistencyCheckDegree(32);
gap> PrimGrpConsistencyCheckDegree(33);
gap> PrimGrpConsistencyCheckDegree(34);
gap> PrimGrpConsistencyCheckDegree(35);
gap> PrimGrpConsistencyCheckDegree(36);
gap> PrimGrpConsistencyCheckDegree(37);
gap> PrimGrpConsistencyCheckDegree(38);
gap> PrimGrpConsistencyCheckDegree(39);
gap> PrimGrpConsistencyCheckDegree(40);
gap> PrimGrpConsistencyCheckDegree(41);
gap> PrimGrpConsistencyCheckDegree(42);
gap> PrimGrpConsistencyCheckDegree(43);
gap> PrimGrpConsistencyCheckDegree(44);
gap> PrimGrpConsistencyCheckDegree(45);
gap> PrimGrpConsistencyCheckDegree(46);
gap> PrimGrpConsistencyCheckDegree(47);
gap> PrimGrpConsistencyCheckDegree(48);
gap> PrimGrpConsistencyCheckDegree(49);
gap> PrimGrpConsistencyCheckDegree(50);
gap> PrimGrpConsistencyCheckDegree(51);
gap> PrimGrpConsistencyCheckDegree(52);
gap> PrimGrpConsistencyCheckDegree(53);
gap> PrimGrpConsistencyCheckDegree(54);
gap> PrimGrpConsistencyCheckDegree(55);
gap> PrimGrpConsistencyCheckDegree(56);
gap> PrimGrpConsistencyCheckDegree(57);
gap> PrimGrpConsistencyCheckDegree(58);
gap> PrimGrpConsistencyCheckDegree(59);
gap> PrimGrpConsistencyCheckDegree(60);
gap> PrimGrpConsistencyCheckDegree(61);
gap> PrimGrpConsistencyCheckDegree(62);
gap> PrimGrpConsistencyCheckDegree(63);
gap> PrimGrpConsistencyCheckDegree(64);
gap> PrimGrpConsistencyCheckDegree(65);
gap> PrimGrpConsistencyCheckDegree(66);
gap> PrimGrpConsistencyCheckDegree(67);
gap> PrimGrpConsistencyCheckDegree(68);
gap> PrimGrpConsistencyCheckDegree(69);
gap> PrimGrpConsistencyCheckDegree(70);
gap> PrimGrpConsistencyCheckDegree(71);
gap> PrimGrpConsistencyCheckDegree(72);
gap> PrimGrpConsistencyCheckDegree(73);
gap> PrimGrpConsistencyCheckDegree(74);
gap> PrimGrpConsistencyCheckDegree(75);
gap> PrimGrpConsistencyCheckDegree(76);
gap> PrimGrpConsistencyCheckDegree(77);
gap> PrimGrpConsistencyCheckDegree(78);
gap> PrimGrpConsistencyCheckDegree(79);
gap> PrimGrpConsistencyCheckDegree(80);
gap> PrimGrpConsistencyCheckDegree(81);
gap> PrimGrpConsistencyCheckDegree(82);
gap> PrimGrpConsistencyCheckDegree(83);
gap> PrimGrpConsistencyCheckDegree(84);
gap> PrimGrpConsistencyCheckDegree(85);
gap> PrimGrpConsistencyCheckDegree(86);
gap> PrimGrpConsistencyCheckDegree(87);
gap> PrimGrpConsistencyCheckDegree(88);
gap> PrimGrpConsistencyCheckDegree(89);
gap> PrimGrpConsistencyCheckDegree(90);
gap> PrimGrpConsistencyCheckDegree(91);
gap> PrimGrpConsistencyCheckDegree(92);
gap> PrimGrpConsistencyCheckDegree(93);
gap> PrimGrpConsistencyCheckDegree(94);
gap> PrimGrpConsistencyCheckDegree(95);
gap> PrimGrpConsistencyCheckDegree(96);
gap> PrimGrpConsistencyCheckDegree(97);
gap> PrimGrpConsistencyCheckDegree(98);
gap> PrimGrpConsistencyCheckDegree(99);

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "consistency99.tst", 1);
