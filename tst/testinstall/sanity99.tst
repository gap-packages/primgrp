gap> START_TEST("sanity.tst");

#
# Disable warnings which depend on Conway Polynomial databases 
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
# Define a function to check the primitive groups of degree n
#
gap> ReadPackage("primgrp", "tst/testutils.g");
true
gap> checkdegree(2);
gap> checkdegree(3);
gap> checkdegree(4);
gap> checkdegree(5);
gap> checkdegree(6);
gap> checkdegree(7);
gap> checkdegree(8);
gap> checkdegree(9);
gap> checkdegree(10);
gap> checkdegree(11);
gap> checkdegree(12);
gap> checkdegree(13);
gap> checkdegree(14);
gap> checkdegree(15);
gap> checkdegree(16);
gap> checkdegree(17);
gap> checkdegree(18);
gap> checkdegree(19);
gap> checkdegree(20);
gap> checkdegree(21);
gap> checkdegree(22);
gap> checkdegree(23);
gap> checkdegree(24);
gap> checkdegree(25);
gap> checkdegree(26);
gap> checkdegree(27);
gap> checkdegree(28);
gap> checkdegree(29);
gap> checkdegree(30);
gap> checkdegree(31);
gap> checkdegree(32);
gap> checkdegree(33);
gap> checkdegree(34);
gap> checkdegree(35);
gap> checkdegree(36);
gap> checkdegree(37);
gap> checkdegree(38);
gap> checkdegree(39);
gap> checkdegree(40);
gap> checkdegree(41);
gap> checkdegree(42);
gap> checkdegree(43);
gap> checkdegree(44);
gap> checkdegree(45);
gap> checkdegree(46);
gap> checkdegree(47);
gap> checkdegree(48);
gap> checkdegree(49);
gap> checkdegree(50);
gap> checkdegree(51);
gap> checkdegree(52);
gap> checkdegree(53);
gap> checkdegree(54);
gap> checkdegree(55);
gap> checkdegree(56);
gap> checkdegree(57);
gap> checkdegree(58);
gap> checkdegree(59);
gap> checkdegree(60);
gap> checkdegree(61);
gap> checkdegree(62);
gap> checkdegree(63);
gap> checkdegree(64);
gap> checkdegree(65);
gap> checkdegree(66);
gap> checkdegree(67);
gap> checkdegree(68);
gap> checkdegree(69);
gap> checkdegree(70);
gap> checkdegree(71);
gap> checkdegree(72);
gap> checkdegree(73);
gap> checkdegree(74);
gap> checkdegree(75);
gap> checkdegree(76);
gap> checkdegree(77);
gap> checkdegree(78);
gap> checkdegree(79);
gap> checkdegree(80);
gap> checkdegree(81);
gap> checkdegree(82);
gap> checkdegree(83);
gap> checkdegree(84);
gap> checkdegree(85);
gap> checkdegree(86);
gap> checkdegree(87);
gap> checkdegree(88);
gap> checkdegree(89);
gap> checkdegree(90);
gap> checkdegree(91);
gap> checkdegree(92);
gap> checkdegree(93);
gap> checkdegree(94);
gap> checkdegree(95);
gap> checkdegree(96);
gap> checkdegree(97);
gap> checkdegree(98);
gap> checkdegree(99);

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "sanity.tst", 1);
