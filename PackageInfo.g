#############################################################################
##  
##  PackageInfo.g for the package `PrimGrp'
##
##  (created from the template file in the `Example' package)
##  

SetPackageInfo( rec(

PackageName := "PrimGrp",
Subtitle := "GAP Primitive Permutation Groups Library",
Version := "4.0.1",
Date := "24/09/2025",
License := "GPL-2.0-or-later",

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
    WWWHome      := "https://www.math.colostate.edu/~hulpke/",
    Place        := "Fort Collins, CO",
    Institution  := Concatenation( [
      "Department of Mathematics, ",
      "Colorado State University",
      ] )
    ),
  rec( 
    LastName      := "Konovalov",
    FirstNames    := "Olexandr",
    IsAuthor      := false,
    IsMaintainer  := true,
    Email         := "obk1@st-andrews.ac.uk",
    WWWHome       := "https://olexandr-konovalov.github.io/",
    PostalAddress := Concatenation( [
                     "School of Computer Science\n",
                     "University of St Andrews\n",
                     "Jack Cole Building, North Haugh,\n",
                     "St Andrews, Fife, KY16 9SX, Scotland" ] ),
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
     ),  

  rec(  # contributed groups of order 4096-8191
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jesse",
    LastName := "Lansdown",
    WWWHome := "https://www.jesselansdown.com",
    Email := "jesse.lansdown@canterbury.ac.nz",
    PostalAddress := Concatenation(
               "Jesse Lansdown\n",
               "School of Mathematics and Statistics\n",
               "University of Canterbury\n",
               "Christchurch 8140\n",
               "New Zealand"),
    Place := "Christchurch",
    Institution := "University of Canterbury",
  ),

  rec(  # contributed groups of order <= 4095
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

  rec(  # contributed groups of order <= 4095
    LastName      := "Russell",
    FirstNames    := "Christopher",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "chriswgrussell@hotmail.co.uk",
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
  permutation groups of degree &lt; 8192.",

PackageDoc := rec(
  BookName  := "primgrp",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "GAP Primitive Permutation Groups Library",
),


Dependencies := rec(
  GAP := "4.10.0",
  NeededOtherPackages := [["GAPDoc", "1.5"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
                      
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := ["primitive permutation group"],

AutoDoc := rec(
    entities := rec(
        VERSION := ~.Version,
        RELEASEDATE := function(date)
          local day, month, year, allMonths;
          day := Int(date{[1,2]});
          month := Int(date{[4,5]});
          year := Int(date{[7..10]});
          allMonths := [ "January", "February", "March", "April", "May", "June", "July",
                         "August", "September", "October", "November", "December"];
          return Concatenation(String(day)," ", allMonths[month], " ", String(year));
        end(~.Date),
        RELEASEYEAR := ~.Date{[7..10]},
    ),
),

));
