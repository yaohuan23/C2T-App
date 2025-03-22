#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
source activate yaohuan23
usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:
  ${NAME}  <input.fastq.gz> 
  record all the fastq.gz file location in a absolute formate 
EOF
}

OUTDIR=$(pwd -P)
if [[ "$1" ]]; then
 if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 1
else
OUTDIR=$1
fi
fi

WRKDIR=$(pwd -P)
errprog=""

SCRIPTARGS="$@"
# main script block
PROMOTER="/public23/home/sca2382/Ref/Annotations/mouse/gencode.vM25.annotation.promoter200.bed"
piprmdup(){
  echo "sample $1 is being treated"
  if [ ! -f $1 ];then
  continue
  fi
  dir=${1%/*}
  filename1=${1%%_*}
  filename=${filename1##*/}
  cd ${dir}
#get the m7G sites by non-template G
#the forward genes
samtools sort -@32 -n  -T ${filename}  $1 -o ${filename}_sorted.bam
samtools fixmate -m ${filename}_sorted.bam ${filename}_sorted-fix.bam
samtools sort -@32 -T ${filename} ${filename}_sorted-fix.bam -o ${filename}_sorted-fix-sorted.bam
samtools markdup -r -@32 ${filename}_sorted-fix-sorted.bam  ${filename}_sorted-fix-sorted-dup.bam
mv ${filename}_sorted-fix-sorted-dup.bam $1
}
pipeline() { 
  echo "sample $1 is being treated"
  if [ ! -f $1 ];then
  continue
  fi
  dir=${1%/*}
  filename1=${1%%_*}
  filename=${filename1##*/}
  cd ${dir}
#get the m7G sites by non-template G
#the forward genes
samtools view -H $1 > ${filename}-Cap-Plus-4G.sam
#samtools view -H $1 > ${filename}-noGend.record-Plus.sam
#samtools view -H $1 > ${filename}-Am.record-Plus.sam
samtools view -f64 -F4 -q 40 $1|grep "XS:A:+"|awk '$6~/^[0-9]+S[1-9][0-9]{1,}M/{print $0}'|awk 'match($6,/^[0-9]+S/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,1,num)~/^[GN]+$/){print $0}}' >> ${filename}-Cap-Plus-4G.sam
#samtools view -f64 -F4 -q 40 $1|grep "XS:A:+"|awk '$6~/^[0-9]+S[1-9][0-9]{1,}M/{print $0}'|awk 'match($6,/^[0-9]+S/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,1,num)~/^[GN]+$/ && substr($10,num+1,1)~/^[ATC]+$/){print $0}}' >> ${filename}-noGend.record-Plus.sam
#samtools view -f64 -F4 -q 40 $1|grep "XS:A:+"|awk '$6~/^[0-9]+S[1-9][0-9]{1,}M/{print $0}'|awk 'match($6,/^[0-9]+S/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,1,num)~/^[GN]+$/ && substr($10,num+1,1)~/^[A]+$/){print $0}}' >> ${filename}-Am.record-Plus.sam

samtools view -bh ${filename}-Cap-Plus-4G.sam > ${filename}-Cap-Plus-4G.bam
#rm ${filename}-Cap-Plus-4G.sam
bedtools bamtobed -i ${filename}-Cap-Plus-4G.bam > ${filename}-Cap-Plus-4G.bed
#rm ${filename}-Cap-Plus-4G.bed
awk '{OFS="\t";$3=$2;print $0}' ${filename}-Cap-Plus-4G.bed > ${filename}-temp
mv ${filename}-temp ${filename}-Cap-Plus-4G.bed
cut -f1-3,6  ${filename}-Cap-Plus-4G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-Plus-m7G.bed
#rm  ${filename}-Cap-Plus-4G.bed
sort -k1,1V ${filename}-Plus-m7G.bed > ${filename}-temp
mv ${filename}-temp ${filename}-Plus-m7G.bed

#samtools view -bh ${filename}-noGend.record-Plus.sam > ${filename}-noGend.record-Plus.bam
#rm ${filename}-noGend.record-Plus.sam
#bedtools bamtobed -i ${filename}-noGend.record-Plus.bam > ${filename}-noGend.record-Plus.bed
#awk '{OFS="\t";$3=$2;print $0}' ${filename}-noGend.record-Plus.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-noGend.record-Plus.bed
#cut -f1-3,6  ${filename}-noGend.record-Plus.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-noGend.record-Plus-m7G.bed
#rm  ${filename}-Cap-Plus-4G.bed
#sort -k1,1V ${filename}-noGend.record-Plus-m7G.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-noGend.record-Plus-m7G.bed


#samtools view -bh ${filename}-Am.record-Plus.sam > ${filename}-Am.record-Plus.bam
#rm ${filename}-noGend.record-Plus.sam
#bedtools bamtobed -i ${filename}-Am.record-Plus.bam > ${filename}-Am.record-Plus.bed
#awk '{OFS="\t";$3=$2;print $0}' ${filename}-Am.record-Plus.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-Am.record-Plus.bed
#cut -f1-3,6  ${filename}-Am.record-Plus.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-Am.record-Plus-m7G.bed
#rm  ${filename}-Cap-Plus-4G.bed
#sort -k1,1V ${filename}-Am.record-Plus-m7G.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-Am.record-Plus-m7G.bed


#The reverse genes
samtools view -H $1 > ${filename}-Cap-minus-4G.sam
samtools view -H $1 > ${filename}-noGend.record-minus.sam
samtools view -H $1 > ${filename}-Am.record-minus.sam
samtools view -f64 -F4 -q 40 $1|grep "XS:A:-"|awk '$6~/[1-9][0-9]{1,}M[0-9]+S$/{print $0}' |awk 'match($6,/[0-9]+S$/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,length($10)-num+1,num)~/^[CN]+$/){print $0}}' >> ${filename}-Cap-minus-4G.sam
samtools view -f64 -F4 -q 40 $1|grep "XS:A:-"|awk '$6~/[1-9][0-9]{1,}M[0-9]+S$/{print $0}' |awk 'match($6,/[0-9]+S$/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,length($10)-num+1,num)~/^[CN]+$/ && substr($10,length($10)-num,1)~/^[ATG]+$/){print $0}}' >> ${filename}-noGend.record-minus.sam
samtools view -f64 -F4 -q 40 $1|grep "XS:A:-"|awk '$6~/[1-9][0-9]{1,}M[0-9]+S$/{print $0}' |awk 'match($6,/[0-9]+S$/){num=substr($6, RSTART, RLENGTH-1);if(substr($10,length($10)-num+1,num)~/^[CN]+$/ && substr($10,length($10)-num,1)~/^[T]+$/){print $0}}' >> ${filename}-Am.record-minus.sam
samtools view -bh ${filename}-Cap-minus-4G.sam > ${filename}-Cap-minus-4G.bam
#rm ${filename}-Cap-minus-4G.sam
bedtools bamtobed -i ${filename}-Cap-minus-4G.bam > ${filename}-Cap-minus-4G.bed
awk '{OFS="\t";$2=$3;print $0}' ${filename}-Cap-minus-4G.bed > ${filename}-temp
mv ${filename}-temp ${filename}-Cap-minus-4G.bed
cut -f1-3,6  ${filename}-Cap-minus-4G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-minus-m7G.bed
#rm ${filename}-Cap-minus-4G.bed
sort -k1,1V ${filename}-minus-m7G.bed > ${filename}-temp
mv ${filename}-temp ${filename}-minus-m7G.bed


#samtools view -bh ${filename}-noGend.record-minus.sam > ${filename}-noGend.record-minus.bam
#rm ${filename}-noGend.record-minus.sam
#bedtools bamtobed -i ${filename}-noGend.record-minus.bam > ${filename}-noGend.record-minus.bed
#awk '{OFS="\t";$2=$3;print $0}' ${filename}-noGend.record-minus.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-noGend.record-minus.bed
#cut -f1-3,6  ${filename}-noGend.record-minus.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-noGend.record-minus-m7G.bed
#rm ${filename}-Cap-minus-4G.bed
#sort -k1,1V ${filename}-noGend.record-minus-m7G.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-noGend.record-minus-m7G.bed


#samtools view -bh ${filename}-Am.record-minus.sam > ${filename}-Am.record-minus.bam
#rm ${filename}-noGend.record-minus.sam
#bedtools bamtobed -i ${filename}-Am.record-minus.bam > ${filename}-Am.record-minus.bed
#awk '{OFS="\t";$2=$3;print $0}' ${filename}-Am.record-minus.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-Am.record-minus.bed
#cut -f1-3,6  ${filename}-Am.record-minus.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > ${filename}-Am.record-minus-m7G.bed
#rm ${filename}-Cap-minus-4G.bed
#sort -k1,1V ${filename}-Am.record-minus-m7G.bed > ${filename}-temp
#mv ${filename}-temp ${filename}-Am.record-minus-m7G.bed
}
pipelineOne() {
while read line;
do
  echo "sample $line is being treated"
  if [ ! -f $line ];then
  continue
  fi
  dir=${line%/*}
  filename1=${line%%_*}
  filename=${filename1##*/}
  normalBam=${line/4G/3G} 
  cd ${dir}
#Generate the Cap-bam from 3G data
samtools view -bh -f64 -F4 -q 40 ${normalBam} > R1-${filename}-3G.bam

#Generate the Cap-bam from 4G data
samtools view -bh -f64 -F4 -q 40 ${line} > R1-${filename}-4G.bam

bedtools bamtobed -i R1-${filename}-3G.bam > R1-${filename}-3G.bed
bedtools bamtobed -i R1-${filename}-4G.bam > R1-${filename}-4G.bed
awk '{OFS="\t"; if($6=="+"){$3=$2;print $0}}' R1-${filename}-3G.bed > Plus-R1-${filename}-3G.bed
awk '{OFS="\t"; if($6=="-"){$2=$3;print $0}}' R1-${filename}-3G.bed > Minus-R1-${filename}-3G.bed

awk '{OFS="\t"; if($6=="+"){$3=$2;print $0}}' R1-${filename}-4G.bed > Plus-R1-${filename}-4G.bed
awk '{OFS="\t"; if($6=="-"){$2=$3;print $0}}' R1-${filename}-4G.bed > Minus-R1-${filename}-4G.bed
#intersectBed -a ${filename}-temp -b ${filename}-Plus-m7G.bed > R1-${filename}-Plus-4G.bed

awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}}' ${filename}-Plus-m7G.bed Plus-R1-${filename}-3G.bed > Plus-R1-${filename}-3G-m7G.bed
awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}}' ${filename}-minus-m7G.bed Minus-R1-${filename}-3G.bed > Minus-R1-${filename}-3G-m7G.bed

cut -f1-3,6  Plus-R1-${filename}-3G-m7G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > All-${filename}-3GPlus-m7G.bedGraph
cut -f1-3,6  Minus-R1-${filename}-3G-m7G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > All-${filename}-3Gminus-m7G.bedGraph

sort -k1,1V All-${filename}-3GPlus-m7G.bedGraph > temp
mv temp All-${filename}-3GPlus-m7G.bedGraph

sort -k1,1V All-${filename}-3Gminus-m7G.bedGraph > temp
mv temp All-${filename}-3Gminus-m7G.bedGraph


awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}}' ${filename}-Plus-m7G.bed Plus-R1-${filename}-4G.bed > Plus-R1-${filename}-4G-m7G.bed
awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}}' ${filename}-minus-m7G.bed Minus-R1-${filename}-4G.bed > Minus-R1-${filename}-4G-m7G.bed

cut -f1-3,6  Plus-R1-${filename}-4G-m7G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > All-${filename}-4GPlus-m7G.bedGraph
cut -f1-3,6  Minus-R1-${filename}-4G-m7G.bed|awk '{$5=$1":"$2;print $0}'|sort -k6|uniq -c |sed 's/^[ ]\{1,\}//g'|sed 's/ /\t/g'|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6,$5}' > All-${filename}-4Gminus-m7G.bedGraph

sort -k1,1V All-${filename}-4GPlus-m7G.bedGraph > temp
mv temp All-${filename}-4GPlus-m7G.bedGraph

sort -k1,1V All-${filename}-4Gminus-m7G.bedGraph > temp
mv temp All-${filename}-4Gminus-m7G.bedGraph

cat All-${filename}-3GPlus-m7G.bedGraph All-${filename}-3Gminus-m7G.bedGraph |awk 'BEGIN{OFS="\t"}{if($6=="-"){$4=$4*-1;print $0}else{print $0}}' > All-${filename}-3GCap.bedGraph
cat All-${filename}-4GPlus-m7G.bedGraph All-${filename}-4Gminus-m7G.bedGraph |awk 'BEGIN{OFS="\t"}{if($6=="-"){$4=$4*-1;print $0}else{print $0}}' > All-${filename}-4GCap.bedGraph
done 
}

pipelineCall() {
  echo "sample $line is being treated"
  if [ ! -f $line ];then
  continue
  fi
  dir=${line%/*}
  filename1=${line%%_*}
  filename=${filename1##*/}
  normalBam=${line/4G/3G} 
  cd ${dir}
   awk 'BEGIN{OFS="\t"}{print $1,$6,$2,$4}' ${filename}-Plus-perfect.bedgraph |sort -k1,1V -k3,3n -k2 > ${filename}-Plus.prePeak
   ~/script/paraclu-9/paraclu 30 ${filename}-Plus.prePeak > ${filename}-Plus.Peak
  sh ~/script/paraclu-9/paraclu-cut.sh ${filename}-Plus.Peak > final-${filename}-Plus.Peak
  rm ${filename}-Plus.prePeak ${filename}-Plus.Peak
  
   awk 'BEGIN{OFS="\t"}{print $1,$6,$2,-1*$4}' ${filename}-Minus-perfect.bedgraph |sort -k1,1V -k3,3n -k2 > ${filename}-minus.prePeak
   ~/script/paraclu-9/paraclu 30 ${filename}-minus.prePeak > ${filename}-minus.Peak
  sh ~/script/paraclu-9/paraclu-cut.sh ${filename}-minus.Peak > final-${filename}-minus.Peak

  cat final-${filename}-Plus.Peak final-${filename}-minus.Peak |awk 'BEGIN{OFS="\t"} $8>0{print $0}'|awk 'BEGIN{OFS="\t";nowChr="";nowStrand="";nowLeft=0;nowRight=0} {if(nowChr==$1 && nowRight+200>$3){nowRight=$4;$3=nowLeft}else{print nowChr,nowLeft,nowRight,nowStrand;nowChr=$1;nowStrand=$2;nowLeft=$3;nowRight=$4} }'|sed '1d' > final-${filename}.Peak
  rm final-${filename}-minus.Peak final-${filename}-Plus.Peak
  rm ${filename}-minus.prePeak ${filename}-minus.Peak 
}

piplineFilter(){
  echo "sample $line is being treated"
  if [ ! -f $line ];then
  continue
  fi
  dir=${line%/*}
  filename1=${line%%_*}
  filename=${filename1##*/}
  normalBam=${line/4G/3G} 
  cd ${dir}
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' ${filename}-noGend.record-minus-m7G.bed  ${filename}-minus-m7G.bed > ${filename}-temp.txt
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' ${filename}-Am.record-minus-m7G.bed  ${filename}-temp.txt > ${filename}-temp12.txt
mv ${filename}-temp12.txt  ${filename}-temp.txt
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' All-${filename}-4Gminus-m7G.bedGraph   ${filename}-temp.txt > ${filename}-temp1.txt
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' All-${filename}-3Gminus-m7G.bedGraph   ${filename}-temp1.txt > ${filename}-temp2.txt
awk 'BEGIN{OFS="\t"} $9/($10+$9) >0.5{print $1,$2,$3,$4*-1,$5,$6}' ${filename}-temp2.txt  > ${filename}-Minus-perfect.bedgraph
#awk 'BEGIN{OFS="\t"} $9/($10+$9) >0.5{print $1,$2,$3,$9*-1,$5,$6}' ${filename}-temp2.txt  > ${filename}-AmMinus-perfect.bedgraph
  

#  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' ${filename}-noGend.record-Plus-m7G.bed  ${filename}-Plus-m7G.bed > ${filename}-temp.txt
#  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' ${filename}-Am.record-Plus-m7G.bed   ${filename}-temp.txt > ${filename}-temp12.txt
mv ${filename}-temp12.txt  ${filename}-temp.txt 
 
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' All-${filename}-4GPlus-m7G.bedGraph   ${filename}-temp.txt > ${filename}-temp1.txt
  awk 'BEGIN{OFS="\t"}NR==FNR{x[$1"-"$2]=$4} NR>FNR{if(x[$1"-"$2]){print $0,x[$1"-"$2]}else{print $0,0}}' All-${filename}-3GPlus-m7G.bedGraph   ${filename}-temp1.txt > ${filename}-temp2.txt
awk 'BEGIN{OFS="\t"} $9/($10+$9) >0.5{print $1,$2,$3,$4,$5,$6}' ${filename}-temp2.txt  > ${filename}-Plus-perfect.bedgraph
#awk 'BEGIN{OFS="\t"} $9/($10+$9) >0.5{print $1,$2,$3,$9,$5,$6}' ${filename}-temp2.txt  > ${filename}-AmPlus-perfect.bedgraph

cat ${filename}-Minus-perfect.bedgraph ${filename}-Plus-perfect.bedgraph > ${filename}-All-perfect.bedgraph
#cat ${filename}-AmMinus-perfect.bedgraph ${filename}-AmPlus-perfect.bedgraph > ${filename}-Am-All-perfect.bedgraph

}

pipenormalize() {
  echo "sample $line is being treated"
  if [ ! -f $line ];then
  continue
  fi
  dir=${line%/*}
  filename1=${line%%_*}
  filename=${filename1##*/}
  normalBam=${line/4G/3G}
  cd ${dir}
  samtools view -c R1-${filename}-3G.bam > ${filename}-counts.txt
  samtools view -c R1-${filename}-4G.bam >> ${filename}-counts.txt
  awk 'BEGIN{a=0}{a+=$1}END{print a}'  ${filename}-counts.txt > ${filename}-countsok.txt
  rm  ${filename}-counts.txt
  echo "awk 'BEGIN{OFS=\"\\t\";getline Counts < \"${filename}-countsok.txt\"}{\$4=\$4*1000*1000000/Counts;print \$0}' ${filename}-All-perfect.bedgraph > CPM-${filename}-All-perfect.bedgraph " > ${filename}-temp.sh
  sh ${filename}-temp.sh

  echo "awk 'BEGIN{OFS=\"\\t\";getline Counts < \"${filename}-countsok.txt\"}{\$4=\$4*1000*1000000/Counts;print \$0}' ${filename}-Am-All-perfect.bedgraph > CPM-${filename}-Am-All-perfect.bedgraph " > ${filename}-temp.sh
  sh ${filename}-temp.sh
}

pipelinotCapBam() {
while read line;
do
  echo "sample $line is being treated"
  if [ ! -f $line ];then
  continue
  fi
  dir=${line%/*}
  filename1=${line%%_*}
  filename=${filename1##*/}
  normalBam=${line/4G/3G}
  cd ${dir}
#Generate the notCap-bam from 3G data
awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}else{print $0>"noGap.bed"}}' ${filename}-Plus-m7G.bed Plus-R1-${filename}-3G.bed > Plus-R1-${filename}-3G-m7G.bed
mv noGap.bed Plus-R1-${filename}-3G-notCap.bed
awk 'BEGIN{OFS="\t"} NR==FNR{x[$1"-"$2]=1}NR>FNR{if(x[$1"-"$2]==1){print $0}else{print $0>"noGap.bed"}}' ${filename}-minus-m7G.bed Minus-R1-${filename}-3G.bed > Minus-R1-${filename}-3G-m7G.bed
mv noGap.bed Minus-R1-${filename}-3G-m7G.bed
cat Minus-R1-${filename}-3G-m7G.bed Plus-R1-${filename}-3G-notCap.bed |cut -f4 > ${filename}-3G-noCap.record
samtools view -H R1-${filename}-3G.bam > R1-${filename}-3G-noCap.sam
samtools view R1-${filename}-3G.bam > R1-${filename}-3G.presam
awk 'NR==FNR{x[$1]=1}NR>FNR{if(x[$1"/1"]){print $0}}' ${filename}-3G-noCap.record R1-${filename}-3G.presam >> R1-${filename}-3G-noCap.sam
samtools view -bh R1-${filename}-3G-noCap.sam > R1-${filename}-3G-noCap.bam
done
}

while read line;
do
pipeline ${line}
done < $1

wait
pipelineOne < $1
wait
pipelinotCapBam < $1
wait

while read line;
do
piplineFilter ${line}
pipenormalize ${line} 
pipelineCall ${line} 
done < $1
