trns=0.005;
snrt=0;

Pr=[1-trns-snrt,trns,snrt;...
    trns,1-2*trns,trns;...
    snrt,trns,1-trns-snrt];