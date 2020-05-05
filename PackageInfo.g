#############################################################################
##  
##  PackageInfo.g for the package `PrimGrp'
##
##  (created from the template file in the `Example' package)
##  

SetPackageInfo( rec(

PackageName := "PrimGrp",
Subtitle := "GAP Primitive Permutation Groups Library",
Version := "3.4.1",
Date := "05/05/2020",
License := "GPL-2.0-or-later",
##  <#GAPDoc Label="PKGVERSIONDATA">
##  <!ENTITY VERSION "3.4.1">
##  <!ENTITY RELEASEDATE "5 May 2020">
##  <!ENTITY RELEASEYEAR "2020">
##  <#/GAPDoc>

PackageWWWHome :=
  Concatenation( "https://gap-packages.github.io/",
      LowercaseString( ~.PackageName ), "/" ),

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", LowercaseString(~.PackageName) ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),

ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", LowercaseString(~.PackageName), "-", ~.Version ),

ArchiveFormats := ".tar.gz",

Persons := [
  rec(
    LastName     := "Hulpke",
    FirstNames   := "Alexander",
    IsAuthor     := true,
    IsMaintainer := true,
    Email        := "hulpke@math.colostate.edu",
    WWWHome      := "http://www.math.colostate.edu/~hulpke",
    Place        := "Fort Collins, CO",
    Institution  := Concatenation( [
      "Department of Mathematics, ",
      "Colorado State University",
      ] )
    ),
  rec( 
    LastName      := "Konovalov",
    FirstNames    := "Alexander",
    IsAuthor      := false,
    IsMaintainer  := true,
    Email         := "alexk@mcs.st-andrews.ac.uk",
    WWWHome       := "http://www.cs.st-andrews.ac.uk/~alexk/",
    PostalAddress := Concatenation( [
                     "School of Computer Science\n",
                     "University of St Andrews\n",
                     "Jack Cole Building, North Haugh,\n",
                     "St Andrews, Fife, KY16 9SX, Scotland" ] ),
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
     ),  
  rec( 
    LastName      := "Roney-Dougal",
    FirstNames    := "Colva M.",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "colva.roney-dougal@st-andrews.ac.uk",
    WWWHome       := "http://www-groups.mcs.st-and.ac.uk/~colva/",
    PostalAddress := Concatenation( [
                     "School of Mathematics and Statistics\n",
                     "University of St Andrews\n",
                     "North Haugh, St Andrews\n",
                     "Fife, KY16 9SS, Scotland" ] ),
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
     ),  
  rec( 
    LastName      := "Russell",
    FirstNames    := "Christopher",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "cr66@st-andrews.ac.uk",
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
     ),  
],

Status := "deposited",

README_URL := 
  Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL := 
  Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),

AbstractHTML := 
  "The <span class=\"pkgname\">PrimGrp</span> package provides the library \
  of primitive permutation groups which includes, up to permutation isomorphism \
  (i.e., up to conjugacy in the corresponding symmetric group), all primitive \
  permutation groups of degree &lt; 4096.",

PackageDoc := rec(
  BookName  := "primgrp",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "GAP Primitive Permutation Groups Library",
),


Dependencies := rec(
  GAP := "4.9.0",
  NeededOtherPackages := [["GAPDoc", "1.5"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
                      
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := ["primitive permutation group"],

));
