library(testthat)
library(pre)

#####
# partykit and earth is loaded as failures of tests may be caused by the version
# of either package

library(earth)
library(partykit)
print(sessionInfo())

test_check("pre")
