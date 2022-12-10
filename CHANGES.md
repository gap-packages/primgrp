# CHANGES to the 'primgrp' GAP package

## 3.4.3 (2022-12-10)

  - Fix inconsistencies in `COHORTS_PRIMITIVE_GROUPS` in degrees 1550 and 1575
  - Fix socle type for natural An/Sn in degree >= 2500
  - Ship data files compressed to save disk space
  - Update some links
  - Minor janitorial changes

## 3.4.2 (2022-05-03)

  - Unify notation for names that involve Mathieu groups to use
    M(n), instead of a mix of M(n) and M_n
  - Fix the name of PrimitiveGroup(81,103) [Reported by Andries Brouwer]
  - Fix the names of PrimitiveGroup(i,j) for (i,j) in
    { (100,14), (121,38), (289,35) }
  - Rename PrimitiveGroup(960,7) from OPlus(8, 2) to O+(8, 2)
    to be consistent with all other groups related to O+ or O-
  - Fix names of PrimitiveGroup(1600,i) for i in [17..19]
  - Fix PrimitiveGroup(1140,2) wrongly named Sym(40) instead of Sym(20)

## 3.4.1 (2020-05-05)

  - Fix name of PrimitiveGroup(625,672);

## 3.4.0 (2019-12-03)

  - Fix the names of PrimitiveGroup(125,17) and PrimitiveGroup(81,125)
  - Improve performance of PrimitiveIdentification
  - Minor janitorial changes

## 3.3.2 (2018-10-27)

  - Use consistent names for the Tits group 2F(4,2)' and for the Ree group
    containing it

## 3.3.1 (2018-02-17)

  - Extend the test suite
  - Minor janitorial changes

## 3.3.0 (2017-12-10)

  - Fix ONanScottType for certain groups of degree 3600
  - Minor janitorial changes

## 3.2.0 (2017-11-11)

  - Document AllPrimitiveGroups and OnePrimitiveGroup
  - Extend the test suite

## 3.1.2 (2017-09-08)

  - Minor janitorial changes

## 3.1.1 (2017-09-01)

  - Add documentation for functions on irreducible solvable matrix groups to the manual

## 3.1.0 (2017-08-31)

  - Initial release as a standalone package
