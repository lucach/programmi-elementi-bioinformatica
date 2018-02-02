#!/bin/bash

dir="tests"
datadir="${dir}/.data"
test -d "$datadir" || mkdir -p "$datadir"
test -f "$datadir"/testo.1.fq.gz || wget -O "$datadir"/testo.1.fq.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz 
test -f "$datadir"/testo.2.fq.gz || wget -O "$datadir"/testo.2.fq.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062635.filt.fastq.gz 

for id in $(seq 2)
do
    text="${datadir}/testo.${id}.fq.gz"
    test -f "$text" || exit 1
    test -f "$dir"/output || mkdir -p "$dir"/output
    for t in ${dir}/input/*
    do
	p=$(cat "$t")
	i=$(basename "$t")
	./karp-rabin --text="$text" --pattern="$p" --rounds=3 --min-quality 20 | grep -oP "\d+" | sort -n > "${dir}/output/$i".kr.${id}
        echo "Done" $id $t
    done
done
diff -uaw --strip-trailing-cr --ignore-all-space "${dir}/output/" "${dir}/ok/"
