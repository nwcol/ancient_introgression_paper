A quick test of support for a late contraction of effective size in the 
Denisova lineage- a feature which appears well-supported due to the shape of the
Denisova D+ decay curve. Also provides an illustration of the bootstrap LRT 
model choice methodology that we are developing.

Nested parameters in this test are T_Dbot and N_Dbot, the time of population
contraction and size following contraction. T_Dbot is a boundary parameter 
(to marginalize the contraction model to the null, set T_Dbot to the sample 
time of Denisova) while N_Dbot is not. The asymptotic null distribution of the
LRT statistic is therefore 1/2 * chi2_1 + 1/2 * ch2_2.

The tests have borderline signifiance (0.05 and 0.07), so we will probably
refrain from implementing this feature in the future unless it appears highly
necessary to obtain a good fit.
