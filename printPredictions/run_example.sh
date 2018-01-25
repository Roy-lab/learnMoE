#!/bin/bash
# Runs a test of learnMoE

#trainM=example_input/ips_trainset0
#testM=example_input/ips_testset0
data=example_input/expression_ipsnames_log
k=3
out=example_output 

mkdir -p $out

echo ./learnMoE -m $data -v3 -o $out -l${k} 
echo gdb ./learnMoE 
echo run -m $data -v3 -o $out -l${k}
