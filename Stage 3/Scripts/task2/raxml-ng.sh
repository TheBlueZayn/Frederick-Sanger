#!/bin/bash

raxml-ng --msa ~/snippydata/core.full.aln --model GTR+G --tree pars{10} --prefix RAxML_result --threads 4 --search
