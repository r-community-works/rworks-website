#
# Implementation of Dyson's algorithm for computing the solution of the
# "Find the dud in M coins in n weighings on a balance scale" problem
# for the special case when M = (3^n - 3)/2.
#
# In particular, this can solve 12 coins in 3 measurements.
#
# This implementation uses signed base-3 representations of numbers,
# in contrast to Dyson, who used standard base-3.


##
## Utility functions on signed trinary
##

#' Convert nonnegative `x` to signed trinary.
#'
#' @param x the nonnegative integer to be converted
#' @param ndigits the number of digits in the trinary representation
#' @returns ndigits-length integer vector `r`, where `r[ndigits]` is the ones digit
#'
to_signed_trinary = function(x, ndigits) {
  maxval = (3^ndigits - 1)/2
  if(x > maxval) stop("x out of range")
  r = numeric(ndigits)
  i = 1
  while(x > 0) {
    d = floor(x/3)
    r[i] = x %% 3
    if(r[i]==2) {
      d = d+1
      r[i] = -1
    }
    x = d
    i = i+1
  }
  rev(r) # put ones digit to the rightmost
}

#' Convert trinary back to decimal
#'
#' @param trin integer vector holding the base-3 representation of a number, ones digit last.
#'
to_decimal = function(trin) {
  ndigits = length(trin)
  w = vapply(1:ndigits, function(i) 3^{ndigits-i}, numeric(1))
  x = sum(w*trin) # dot product of w and trinary representation
  x
}

#' Determine if a signed trinary number rotates "forward".
#'
#' A number goes "forward" if the first nonzero difference in digits is 1 mod 3.
#' This progression is -1 -> 0 -> 1 -> -1 ... .
#' Note that (-1) - 1 = -2, and -2 %% 3 = 1
#'
#' @param ptrinary integer vector holding the base-3 representation of a number, ones digit last.
#' @returns TRUE if ptrinary rotates forward, else FALSE. Will err out on a constant vector.
#'
is_forward = function (ptrinary) {
  delta = diff(ptrinary)
  # get only the nonzero diffs
  delta = delta[!(delta==0)]

  if(length(delta)==0)
    stop("constant vector")

  # check if we are rotating forward, and return
  delta[1] %% 3 == 1
}

##
## Supporting functions for the algorithm
##

#'  Get the forward representation of the maximum possible number of coins, to `ndigit` digits.
#'
#' Gets the forward representations of the integers 1:M, where `M = (3^ndigits - 3)/2`. If the
#' representation of the integer `k` does not rotate forward, then the representation
#' of `-k` does. We'll say that if a representation does not rotate forward, then it
#' rotates.backwards.
#'
#' @param ndigits the number of digits in the trinary representation
#' @returns a matrix `Fmat` of dimensions `ncoins` x `ndigits` such that `Fmat[i, ]` represents `i`
#'
get_forward_representation = function(ndigits) {
  ncoins = (3^ndigits - 3)/2
  Fmat = matrix(0, nrow=ncoins, ncol=ndigits)
  for(i in 1:ncoins) {
    trin = to_signed_trinary(i, ndigits)
    if(is_forward(trin))
      Fmat[i, ] = trin
    else
      Fmat[i, ] = -1*trin
  }
  Fmat
}



#' Calculate the positions of each coin for each weighing
#'
#' The forward representation of a coin dictates where it goes for the each weighing.
#' -1 designates the left pan, 1 designates the right pan, and 0 is the table.
#' The ith digit positions the coin for the ith weighing, so the trinary representation
#' should have as many digits as there are weighings (rounds).
#'
#' The function returns a list `C_all` of length the number of rounds.
#' `C_all[i]]` is a data frame with columns "left", "table" and "right". Each column is
#' a vector of coins (integers) to be placed in the appropriate location on the ith weighing.
#'
#' We can precompile this information, as it never varies for a given number of rounds.
#'
#' @param Fmat the matrix of forward representations, as given by `get_forward_representation`.
#' @returns the list `C_all`
#'
compileC = function(Fmat) {
  nrounds = ncol(Fmat)
  ncoins = nrow(Fmat)

  C_all = list()

  for(iround in 1:nrounds) {
    # get the ith digit of all the numbers
    digits = as.numeric(Fmat[ , iround ])
    left = which(digits==-1)
    table = which(digits==0)
    right = which(digits==1)

    C_all[[iround]] = data.frame(left=left, table=table, right=right)
  }

  C_all

}

#' "Public" function to precompile what's needed for solving the `ncoins` in `nrounds` problem.
#'
#' @param ncoins The number of coins to weigh
#' @param nrounds The number of weighings to do
#' @returns a list `Cset` that describes the coin weighing schedule.
#'
precompile = function(ncoins, nrounds) {
  if(ncoins < 3) stop("ncoins must be at least 3")
  M_max = (3^nrounds -3)/2
  if(ncoins > M_max) stop("You will need more rounds for this many coins.")
  if(ncoins!=M_max) stop("I've only implemented ncoins==M_max; sorry.")

  # encode the coins
  Fmat = get_forward_representation(nrounds)
  # calculate the coin partitioning schedule
  Cset = compileC(Fmat)

  Cset
}



#' Weigh a vector of coins according to the weighing schedule
#'
#' `coins` is a vector of coin weights, where at most one coin weighs
#' differently than the others (the "dud"). The outcome of each weighing
#' is the pan that was heavier (if any):
#'
#' * If the left pan is heavier, return -1
#' * If the right pan is heavier, return 1
#' * If the scale balances, return 0
#'
#' The vector `a` of outcomes is the signed trinary representation of the dud.
#' If the representation rotates forward, the dud is heavier than the good coins.
#' If the representation rotates backward, the dud is lighter than the good coins.
#' All zeros means there is no dud.
#'
#' @param coins the vector of coin weights
#' @param Cset the precompiled weighing schedule
#' @returns a vector `a` of outcomes. the same length as `Cset`
#'
do_weighings = function(coins, Cset) {
  nrounds = length(Cset)
  ncoins = length(coins)
  a = numeric(nrounds)

  for(iround in 1:nrounds) {
    Ci = Cset[[iround]]
    leftwt = sum(coins[Ci$left])
    rightwt = sum(coins[Ci$right])
    a[iround] = ifelse(leftwt==rightwt, 0,
                       ifelse(leftwt > rightwt, -1, 1))
  }
  a
}


#' Given a vector of coins, and a number of rounds, find the dud if it exists, and its weight.
#'
#' This is the "public" call to find the solution for a specific vector of coin weights,
#' assuming that the necessary weighing schedule has been precompiled.
#'
#' The result `A` is derived from the vector of weighing outcomes `a`, which
#' is the signed trinary representation of the dud. The absolute value of
#' `to_decimal(a)` gives the index of the dud coin (0 if there is no dud).
#' We define `A` to be positive if the dud is heavier than the good coins,
#' and negative if the dud is lighter than the good coins.
#'
#' Note that `sign(A)` is NOT the same as `sign(to_decimal(a))`, which merely notes
#' whether the trinary representation of a given positive integer rotates forward or backward.
#'
#' @param coins the vector of coin weights
#' @param Cset the precompiled weighing schedule
#' @returns an integer whose magnitude gives the dud location, and whose sign designates the relative weight.
#' `A=0` means there is no dud.
#'
find_dud = function(coins, Cset) {
  ncoins = length(coins)
  precoins =  (3^nrounds -3)/2
  if(ncoins != precoins) stop("Wrong number of coins")

  a = do_weighings(coins, Cset)
  A = to_decimal(a)

  if(A==0) {
    print("No dud.")
    return(A)
  }

  direction = ifelse(is_forward(a) > 0, 1, -1)
  abs(A)*direction
}



