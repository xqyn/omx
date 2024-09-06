library(devtools)
install_github("dasonk/docstring")
library(docstring)



square <- function(x){
    #' Square a number
    #'
    #' Calculates the square of the input
    #'
    #' @param x the input to be squared
    #' @param k the input to be squared
    #' @param f the input to be squared

    return(x^2)
}
docstring(square)




add_numbers <- function(x, y) {

    #' Add two numbers
    #'
    #' This function adds two numbers together.
    #' @param x A numeric value.
    #' @param y A numeric value.
    #' @return The sum of \code{x} and \code{y}.
    #' @examples
    #' add_numbers(3, 5)
    #' @export

  x + y
}

docstring(add_numbers)

library(roxygen2)
roxygenize()