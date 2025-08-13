Labels for features:
a: S->ND pulse
b: S->D pulse
c: MH->N pulse 
d: D-N migration

Labels for models:
0:  null
1:  a       s->nd
2:  b       s->d
3:  c       mh->n
4:  d       d<->n
5:  ab      s->nd, s->d
6:  ac      s->nd, mh->n
7:  ad      s->nd, d<->n
8:  bc      s->d, mh->n
9:  bd      s->d, d<->n
10: cd      mh->n, d<->n
11: abc     s->nd, s->d, mh->n
12: abd     s->nd, s->d, d<->n
13: acd     s->nd, mh->n, d<->n
14: bcd     s->d, mh->n, d<->n
15: abcd    s->nd, s->d, mh->n, d<->n

Here are model labels expressed as a Python dictionary:
model_labels = {
    0: "null",
    1: "S->ND",
    2: "S->D",
    3: "H->N",
    4: "D<->N",
    5: "S->ND,S->D",
    6: "S->ND,H->N",
    7: "S->ND,D<->N",
    8: "S->D,H->N",
    9: "S->D,D<->N",
    10: "H->N,D<->N",
    11: "S->ND,S->D,H->N",
    12: "S->ND,S->D,D<->N",
    13: "S->ND,H->N,D<->N",
    14: "S->D,H->N,D<->N",
    15: "S->ND,S->D,H->N,D<->N"
}