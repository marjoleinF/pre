library(testthat)
library(pre)

#####
# partykit and earth is loaded as failures of tests may be caused by the version
# of either package. Thus, we print the sessionInfo

library(earth)
library(partykit)
print(sessionInfo())

test_check("pre")
