#!/bin/bash
quantile=("0.95" "0.99" "0.999")
linbash=("TETDIP" "PANDIN" "PANSEC" "PANBAL" "PANWCA" "DINSEC" "DINBAL" "DINWCA" "SECWCA" "BALSEC" "BALWCA" "CROPAN" "CRODIN" "CROSEC" "CROBAL" "CROWCA")
    for quant in "${quantile[@]}"; do
    head -n1 data/FstStatsPerContrast.$quant.${linbash[0]}.txt > results/FstStatsPerContrast.$quant.txt
    cut -f1 data/FstHighperGene.Allgenes.$quant.${linbash[0]}.txt > FstHighperGene.Allgenes.$quant.txt
        for lin in "${linbash[@]}"; do
        tail -n1 data/FstStatsPerContrast.$quant.$lin.txt >> results/FstStatsPerContrast.$quant.txt
        tail -n +1 data/FstHighperGene.Allgenes.$quant.$lin.txt | cut -f2 > aa
        paste FstHighperGene.Allgenes.$quant.txt aa > aaa
        mv aaa FstHighperGene.Allgenes.$quant.txt
        done
    mv FstHighperGene.Allgenes.$quant.txt results/FstHighperGene.Allgenes.$quant.txt 
    done

rm aa
