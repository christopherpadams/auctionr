## Resubmission
This is a resubmission. In this version I have:

* Removed trailing spaces from the DESCRIPTION file.

* Added a reference to the methods paper the package is based on to the
  Description field.

* Removed unnecessary console logging from the parameter estimation routine.

NOTE: The first two expressions in the auction_model example run in < 5 seconds
and are sufficient to test the function. The remaining expressions are wrapped
in `\dontrun` because they take a long time to run. These remaining expressions
are important to demonstrate real-world use, and cannot be simplified without
losing their value.
