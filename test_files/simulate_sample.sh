#!/bin/bash

# 设置参数
TOTAL_READS=1000000  # 总共生成1百万条reads
READS_90=$(echo "$TOTAL_READS * 0.9" | bc | cut -d. -f1)
READS_10=$(echo "$TOTAL_READS * 0.1" | bc | cut -d. -f1)

# 组合1: EC590(90%) + K12(10%)
echo "模拟组合1: EC590(90%) + K12(10%)"
wgsim -N $READS_90 e.coli-EC590.fasta tmp_ec590_1.fq tmp_ec590_2.fq
wgsim -N $READS_10 e.coli-K12.fasta tmp_k12_1.fq tmp_k12_2.fq
cat tmp_ec590_1.fq tmp_k12_1.fq > EC590_K12_1.fq
cat tmp_ec590_2.fq tmp_k12_2.fq > EC590_K12_2.fq
rm tmp_*.fq

# 组合2: EC590(90%) + O157(10%)
echo "模拟组合2: EC590(90%) + O157(10%)"
wgsim -N $READS_90 e.coli-EC590.fasta tmp_ec590_1.fq tmp_ec590_2.fq
wgsim -N $READS_10 e.coli-o157.fasta tmp_o157_1.fq tmp_o157_2.fq
cat tmp_ec590_1.fq tmp_o157_1.fq > EC590_O157_1.fq
cat tmp_ec590_2.fq tmp_o157_2.fq > EC590_O157_2.fq
rm tmp_*.fq

# 组合3: K12(90%) + O157(10%)
echo "模拟组合3: K12(90%) + O157(10%)"
wgsim -N $READS_90 e.coli-K12.fasta tmp_k12_1.fq tmp_k12_2.fq
wgsim -N $READS_10 e.coli-o157.fasta tmp_o157_1.fq tmp_o157_2.fq
cat tmp_k12_1.fq tmp_o157_1.fq > K12_O157_1.fq
cat tmp_k12_2.fq tmp_o157_2.fq > K12_O157_2.fq
rm tmp_*.fq

echo "所有组合模拟完成"
