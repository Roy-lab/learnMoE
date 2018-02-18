#!/bin/bash
# Runs a test of learnMoE with tab reading capabilities

#trainM=example_input/ips_trainset0
#testM=example_input/ips_testset0
data=example_input/ips_trainset0.txt
k=3
out=example_output 

mkdir -p $out

echo ./learnMoE -m $data -o $out -l${k} 
echo gdb ./learnMoE 
echo run -m $data -o $out -l${k}
