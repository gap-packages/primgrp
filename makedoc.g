##  this creates the documentation, needs: GAPDoc and AutoDoc packages, pdflatex
##
##  Call this with GAP from within the package directory.

if fail = LoadPackage("AutoDoc", ">= 2025.12.19") then
    Error("AutoDoc 2025.12.19 or newer is required");
fi;

AutoDoc(rec( scaffold := rec( includes := [ "prim.xml", "irredsol.xml" ],
                              bib := "manualbib.xml" ),
             extract_examples := true ));
Exec("mv tst/primgrp01.tst tst/manualexamples/");
