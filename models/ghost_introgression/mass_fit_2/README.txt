Here we fit a reduced set of models.

Model layout.
    S pulses           none    S->ND   S->D    both
H->N pulse?     No     0       1       2       3
                Yes    4       5       6       7


observed ll.
                        S->X pulses?
                        none    S->ND   S->D    S->ND and S->D
H->N pulses?    No      -605    -603    -408    -408
                Yes     -303    -300    -298    -297




Compare to older results:
 model           description          ll  best_fit
     4                 D<->N -604.906812         2
     7           S->ND,D<->N -602.693604        23
     9            S->D,D<->N -408.901920         9
    12      S->ND,S->D,D<->N -408.991943         3

    10            H->N,D<->N -394.134244         4
    13      S->ND,H->N,D<->N -392.461686        10
    14       S->D,H->N,D<->N -386.950477        27
    15 S->ND,S->D,H->N,D<->N -385.428760        27