#
# Implementation of Dyson's algorithm for computing the solution of the
# "Find the dud in M coins in n weighings on a balance scale" problem
# for the general case 3 <= M <= (3^n - 3)/2.
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

#' Increase all digits of a signed trinary number by 1 mod 3
#'
#' @param a number in signed trinary
#' @returns the shifted number
cyclic_shift = function(trinary) {
  shifted = (trinary + 1)
  shifted = ifelse(shifted==2, -1, shifted)
  shifted
}

##
## Supporting functions for the algorithm
##

#' Get the forward representation of the maximum possible number of coins, to `ndigit` digits.
#'
#' Gets the forward representations of the integers 1:M, where `M = (3^ndigits - 3)/2`.
#' If the representation of the integer `k` does not rotate forward, then the representation
#' of `-k` does. We'll say that if a representation does not rotate forward, then it
#' rotates.backwards.
#'
#' @param ndigits the number of digits in the trinary representation
#' @returns a matrix `Fmat` of dimensions `M` x `ndigits` such that `Fmat[i, ]` represents `i`
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


#' Create groups of "weighing triplets" for nrounds weighings
#'
#' Divides the `M_max = (3^nrounds -3)/2` possible labels into groups
#' of three, such that every group represents three coins weighed
#' on the scale. This means that no two numbers in a triplet
#' have the same value digit in any position.
#'
#' @param nrounds the number of weights (and the number of digits)
#' @returns a vector `g = 1:M_max` where `names(g)` gives the group for each integer
#'
create_cyclic_groups = function(nrounds) {
  maxM = (3^nrounds -3)/2
  ngroups = maxM/3
  groups = numeric(maxM)

  for(ig in 1:ngroups) {
    # get the first unassigned id
    ix = which(groups==0)[1]
    for(ishift in 1:3) {
      groups[abs(ix)] = ig
      ixtrin = to_signed_trinary(ix, nrounds)
      ixtrin = cyclic_shift(ixtrin)
      ix = to_decimal(ixtrin)
    }
  }
  groups
}

#' Assign a label to every coin
#'
#' @param ncoins the number of coins
#' @param the number of weighings (and digits)
#' @param the group each label belongs to, as produced by `create_cyclic_groups`
#' @returns a vector `l = 1:ncoins`, where `names(l)` gives the label for each coin
#'
assign_labels = function(ncoins, ndigits, groups) {
  labels = numeric(ncoins)

  nwhole = floor(ncoins/3)
  rem = ncoins %% 3

  for(ng in 1:nwhole) {
    starti = 3*(ng-1)+1
    endi = 3*ng
    labels[starti:endi] = which(groups==ng)
  }

  if(rem > 0) {
    lastgp = which(groups==nwhole+1)
    ix = 3*nwhole + 1
    # not the most efficient way to do this, but it's only 3
    for(x in lastgp) {
      itrin = to_signed_trinary(x, ndigits)
      if((rem==1) && (itrin[1]==0)) {
        labels[ix] = x
      }
      if((rem==2) & (itrin[1] != 0)) {
        labels[ix] = x
        ix = ix+1
      }
    }
  }

  labels

}



#' Calculate the positions of each possible label for each weighing
#'
#' The forward representation of a label dictates where it goes for the each weighing.
#' -1 designates the left pan, 1 designates the right pan, and 0 is the table.
#' The ith digit positions the labeled coin for the ith weighing, so the trinary representation
#' should have as many digits as there are weighings (rounds).
#'
#' The function returns a list `C_all` of length the number of rounds.
#' `C_all[i]]` is a data frame with columns "left", "table" and "right". Each column is
#' a vector of labels (integers) to be placed in the appropriate location on the ith weighing.
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
#' Returns the weighing schedule as calculated by `compileC`, and the coin groups as
#' calculated by `compile_cyclic_groups`
#' @param nrounds The number of weighings to do
#' @returns a list L where L$Cset is the weighing schedule, and L$groups is the group assigment vector
#'
precompile = function(nrounds) {
  # encode the coins
  Fmat = get_forward_representation(nrounds)
  # calculate the coin partitioning schedule
  Cset = compileC(Fmat)

  # get the cyclic groupss for when M <> M_max
  groups = create_cyclic_groups(nrounds)

  list(Cset = Cset, groups = groups)
}


#' From a set of possible labels, get the labels that are actually used
#'
#' The coins are all the coins that we are evaluating. The labels are
#' all the possible labels for a particular set (say, the left pan).
#'
#' If there are fewer than `(3^nrounds -3)/2` coins, then some of the labels
#' in `all_labels` are not used. We want to intersect the labels of our coins
#' with the set of possible labels, and return the resulting set of coins.
#'
#' @param coins a set of labelled coins
#' @param all_labels a set of valid labels to pick from
#' @returns exactly the set of coins that have labels from `all_labels`
#'
get_coins = function(coins, all_labels) {
  validset = intersect(all_labels, names(coins))
  coins[as.character(validset)]
}


#' From a set of coins, pick one that is known to be good.
#' @param coinset a set of labelled coins
#' @param is_good a logical logical vector that is TRUE for coins we know to be good
#' @returns a good coin if one exists, otherwise NULL
#'
find_good_coin = function(coinset, is_good) {
  mystats = is_good[names(coinset)]
  g_ix = which(mystats)

  if(FALSE) {
    print("finding a good coin from")
    print(coinset)
    print("stats")
    print(mystats)
    print("good coins")
    print(g_ix)
    if(length(g_ix) > 0) print(paste("returning", coinset[g_ix[1]] ))
  }

  if(length(g_ix)==0) return(NULL)

  coinset[g_ix[1]]

}

LEFT = 1
RIGHT = 2

# do the nrounds weighing of the pennies
# and return the "signature" of the dud

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
#' @param labels the labels for the coins
#' @returns a vector `a` of outcomes. the same length as `Cset`
#'
do_weighings = function(coins, Cset, labels) {
  nrounds = length(Cset)
  ncoins = length(coins)
  a = numeric(nrounds)
  scale = list()

  names(coins) = as.character(labels)
  is_good = logical(ncoins)  # initialized to FALSE
  names(is_good) = as.character(labels)


  for(iround in 1:nrounds) {
    Ci = Cset[[iround]]

    leftset = get_coins(coins, Ci$left)
    scale[[LEFT]] = leftset
    leftlen = length(scale[[LEFT]])

    rightset = get_coins(coins, Ci$right)
    scale[[RIGHT]] = rightset
    rightlen = length(scale[[RIGHT]])

    tablecoins = get_coins(coins, Ci$table)

    if(FALSE) {
      print("beginning of round")
      print(scale)
      print(tablecoins)
      print(is_good)
    }

    if(leftlen != rightlen) {
      if(leftlen < rightlen) {
        shorter = LEFT
        longer = RIGHT
      } else {
        shorter = RIGHT
        longer = LEFT
      }

      # get a known good coin from the table and
      # add it to the side that's missing a coin
      good_coin = find_good_coin(tablecoins, is_good)
      if(!is.null(good_coin)) {
        scale[[shorter]] = c(scale[[shorter]], good_coin)
      } else {
        # otherwise, move a coin from the larger set to the table
        # which we can implement by just removing a good coin
        # from the longer set.
        good_coin = find_good_coin(scale[[longer]], is_good)
        remaining = setdiff(names(scale[[longer]]), names(good_coin))
        scale[[longer]] = scale[[longer]][remaining]
      }
    }

    if(FALSE){
      print("before weighing")
      print(scale)
    }

    leftwt = sum(scale[[LEFT]])
    rightwt = sum(scale[[RIGHT]])
    a[iround] = ifelse(leftwt==rightwt, 0,
                       ifelse(leftwt > rightwt, -1, 1))
    if(FALSE) {
      print(paste("result of weighing", a[iround]))
    }


    # update the statuses
    if(a[iround]==0) { # the coins on the scale are good
      scalenames = c(names(scale[[LEFT]]), names(scale[[RIGHT]]))
      is_good[scalenames] = TRUE
    } else { # the coins on the table are good
      is_good[names(tablecoins)] = TRUE
    }
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
find_dud = function(coins, precompiled) {
  Cset = precompiled$Cset
  groups = precompiled$groups

  nrounds = length(Cset)
  M_max = (3^nrounds -3)/2
  ncoins = length(coins)

  if(ncoins < 3) stop("Must have at least 3 coins")
  if(ncoins > M_max) stop("You will need more rounds for this many coins.")

  if(ncoins == M_max) { # just run the basic algo
    labels = 1:ncoins
  } else {
    labels = assign_labels(ncoins, nrounds, groups)
  }

  a = do_weighings(coins, Cset, labels)
  A = to_decimal(a)

  if(A==0) {
    print("No dud.")
    return(A)
  }

  badcoin = which(labels==abs(A))

  direction = ifelse(is_forward(a) > 0, 1, -1)
  badcoin*direction
}
