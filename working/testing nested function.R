inside_print <- function(mystring){paste(mystring, "2")}

print_test <- function(string) {

  inside_print <- function(mystring){ paste(mystring, "2")}

  paste("2", inside_print(string))
}

print_test("something")
