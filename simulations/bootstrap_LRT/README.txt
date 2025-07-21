This group of simulations is intended to clarify the correct form of adjustment
for the composite LRT by comparing several adjustment setups to model rejection
rates of fits to bootstrap replicates. These pairs of null/alternate models are 
treated (numbers in parenthesis are the numbers of on-/off-boundaray nested
parameters):

isolation/symmetric migration                           (1, 0)  
symmetric size/asymmetric size (under isolation)        (0, 1)
symmetric migration/asymmetric migration                (0, 1)
equilibrium/contraction                                 (1, 1)
    fixed change time                                   (0, 1)
isolation/admixture pulse                               (1, 1)
    fixed pulse time                                    (1, 0)

Throughout I use these labels for the LRT adjustments:
    0A: evaluation at null MLE marginalized to null parameterization
    0B: evaluation at null MLE with nested parameters set to model 1 MLE
    1: evaluation at model 1 MLE






Cluster IDs:
asymmetric size: 4416112
asymmetric mig: 4416116 4416205
migration: 4416132
contraction: 4417643
pulse: 4417645