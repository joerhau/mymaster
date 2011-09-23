#!/bin/bash
#make -f Makefile.JOERG.PTHREADS.gcc
./raxmlLight-PTHREADS -T 4 -t workin/RAxML_parsimonyTree.START -s workin/prot.phy -m PROTGAMMAWAG -q workin/partitionFile.txt -n Q > out.txt
