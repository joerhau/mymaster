#!/bin/bash
raxmlLight-PTHREADS -T 4 -t RAxML_parsimonyTree.START -s $1.phy -m PROTGAMMAWAG -q $1.part -n Q1 -l

