#############################################################################
##  
##  PackageInfo.g for the package `PrimGrp'
##
##  (created from the template file in the `Example' package)
##  

SetPackageInfo( rec(

PackageName := "PrimGrp",
Subtitle := "GAP Primitive Permutation Groups Library",
Version := "4.1.0-DEV",
Date := "28/07/2026",
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
    PostalAddress := Concatenation( [
                     "Department of Mathematics\n",
                     "Colorado State University\n",
                     "Fort Collins, CO, 80523-1874, USA" ] ),
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

  rec(  # reworked how the groups are stored
    LastName       := "Horn",
    FirstNames     := "Max",
    IsAuthor       := false,
    IsMaintainer   := true,
    Email          := "mhorn@rptu.de",
    WWWHome        := "https://www.quendi.de/math",
    GitHubUsername := "fingolfin",
    PostalAddress  := Concatenation( [
                      "Fachbereich Mathematik\n",
                      "RPTU Kaiserslautern-Landau\n",
                      "Gottlieb-Daimler-Straße 48\n",
                      "67663 Kaiserslautern\n",
                      "Germany" ] ),
    Place          := "Kaiserslautern, Germany",
    Institution    := "RPTU Kaiserslautern-Landau"
     ),

  rec(  # contributed groups of order 4096-8191
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jesse",
    LastName := "Lansdown",
    WWWHome := "https://www.jesselansdown.com",
    Email := "jesse.lansdown@canterbury.ac.nz",
    PostalAddress := Concatenation(
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
  GAP := "4.12.0",
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []
                      
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := ["primitive permutation group"],

AutoDoc := rec(
    TitlePage := rec(
        Abstract := """
            <Index Key="PrimGrp package">&primgrp; package</Index>
            The &GAP; package &primgrp; provides the library of primitive
            permutation groups which includes, up to permutation isomorphism
            (i.e., up to conjugacy in the corresponding symmetric group),
            all primitive permutation groups of degree &lt; 8192.""",
        Copyright := """
            &primgrp; is free software; you can redistribute it and/or modify it
            under the terms of the GNU General Public License as published by
            the Free Software Foundation; either version 2 of the License, or
            (at your option) any later version. For details, see the FSF's own site
            <URL>https://www.gnu.org/licenses/gpl.html</URL>.
            <P/>

            If you publish a result which was partially obtained with the usage of
            &primgrp;, please cite it in the following form:
            <P/>

            A. Hulpke, J. Lansdown, C. Roney-Dougal, C. Russell.
            <E>PrimGrp --- GAP Primitive Permutation Groups Library,
            Version &VERSION;;</E> &RELEASEYEAR;
            (<URL>https://gap-packages.github.io/primgrp/</URL>).
            <P/>

            For where the data itself comes from, and who computed which part
            of it, see Section <Ref Sect="Overview"/>.""",
        Acknowledgements := """
            The conversion of the &GAP; database of primitive permutation groups to
            a separate &GAP; package has been supported by the EPSRC Collaborative
            Computational Project EP/M022641/1 CoDiMa (CCP in the area of Computational
            Discrete Mathematics), <URL>https://www.codima.ac.uk/</URL>.""",
    ),
    entities := rec(
        primgrp := "<Package>PrimGrp</Package>",
    ),
),

));
