
FPA – the four-point aligner
----------------------------

**Compilation**

    $ git clone https://github.com/ariloytynoja/fpa.git
    $ cd fpa/
    $ make

**Usage**

FPA has two use modes, scan and visualisation. For filtering, awk is used.

First, mutation clusters are identified and the best explanation involving a template switch is computed. The results are outputted as a table that is redirected to a file.

    $ ./fpa --scan --pair homo_sapiens.12.74743744.74973891.fas > homo_sapiens.12.74743744.74973891.csv

Second, candidate events are filtered using awk.

    $ awk -F, '($8-$9>12 && $11>=0.95 && $13>=0.95 && $17==0 && $20>5) {print $0}' homo_sapiens.12.74743744.74973891.csv  > one_hit.csv

More complex criteria can be used for filtering. The fields used here
are:

8: sp2_ref     (switch point 2)
9: sp3_ref     (switch point 3)
11: iden_up    (identity upstream)
13: ident_down (identity downstream)
17: masked     (sequence masking: 0=none)
20: sum_mis    (excess mismatches in forward alignments)

Third, interesting cases are visualised.

    $ ./fpa --pair homo_sapiens.12.74743744.74973891.fas --print-file one_hit.csv 

    chr12:74744825-74744838

    Switch process:
    F1: L AATTACACAATTGTTGATAATATTGTGCTGGCCTTGTTCC 1
    F3:                                                          4 TGTTGCATGTCATATCATTACCCTTAACTATGTAATCT R
    RF:   AATTACACAATTGTTGATAATATTGTGCTGGCCTTGTTCCCCATTTAGCTACCTTCATGTTGCATGTCATATCATTACCCTTAACTATGTAATCT
    RR:   TTAATGTGTTAACAACTATTATAACACGACCGGAACAAGGGGTAAATCGATGGAAGTACAACGTACAGTATAGTAATGGGAATTGATACATTAGA
    F2:                                   3 ACAAGGGGTAAATC 2

    Template-switch alignment:
     TATTGTGCTGGCCTTGTTCC|CTAAATGGGGAACA|TGTTGCATGTCATATCATTA
     TATTGTGCTGGCCTTGTTCC|CTAAATGGGGAACA|TGTTGCATGTCATATCATTA

    EPO alignment:
     TATTGTGCTGGCCTTGTTCCCtAaaTgGggAa---CATGTTGCATGTCATATCATTA
     TATTGTGCTGGCCTTGTTCCCCATTTAGCTACCTTCATGTTGCATGTCATATCATTA

