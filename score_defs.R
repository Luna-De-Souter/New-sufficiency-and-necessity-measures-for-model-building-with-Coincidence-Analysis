
# List of score definitions to be benchmarked Version 1
# #####################################################

score_defs <- list(
# control trials, regular consistency, regular coverage:
wRatio_control = ccDef_ratio(),
 
# PA-consistency, AA-coverage:
wRatio_def1 = ccDef_ratio(list(
       cbind(c(1, 0, 0, 0, 0, 0),
       c(1, 0, -99-1i, 0, -99-1i, 0)), #-99-1i => (sy/SY)^-1 = sY/sy
       cbind(c(1, 0, 0, 0, 0, 0),
       c(1, 0, 0, -98-1i, -98-1i, 0)))),

# AAC-consistency, AA-coverage:
wRatio_def2 = ccDef_ratio(list(
       cbind(c(0, 1, 0, 0, 0, 0),
       c(0, 1, -98, 0, -98, 0)),
       cbind(c(1, 0, 0, 0, 0, 0),
       c(1, 0, 0, -98-1i, -98-1i, 0)))),

# PA-consistency, PAC-coverage:
wRatio_def3 = ccDef_ratio(list(
       cbind(c(1, 0, 0, 0, 0, 0),
       c(1, 0, -99-1i, 0, -99-1i, 0)),
       cbind(c(0, 1, 0, 0, 0, 0),
       c(0, 1, 0, -99, -99, 0)))),
    
# AAC-consistency, PAC-coverage:
wRatio_def4 = ccDef_ratio(list(
       cbind(c(0, 1, 0, 0, 0, 0),
       c(0, 1, -98, 0, -98, 0)),
       cbind(c(0, 1, 0, 0, 0, 0),
       c(0, 1, 0, -99, -99, 0))))
)


