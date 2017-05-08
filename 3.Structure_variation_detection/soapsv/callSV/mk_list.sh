#argv1 input all result file
#argv2 file fix

awk '($24>=70||$24==50)&&$15<=30&&$18<=10&&($16/(100-$15))<0.1&&($17/(100-$15))<0.1' $1 > ${1}_confident
awk '!(($24>=70||$24==50)&&$15<=30&&$18<=10&&($16/(100-$15))<0.1&&($17/(100-$15))<0.1)' $1 > ${1}_potential
awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$30}' ${1}_confident > BGI_${2}_confident_indels
awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$30}' ${1}_potential > BGI_${2}_potential_indels
awk '$1=="Deletion"' BGI_${2}_confident_indels > BGI_${2}_confident_indels_deletion
awk '$1=="Deletion"' BGI_${2}_potential_indels > BGI_${2}_potential_indels_deletion
awk '$1=="Insertion"' BGI_${2}_confident_indels > BGI_${2}_confident_indels_insertion
awk '$1=="Insertion"' BGI_${2}_potential_indels > BGI_${2}_potential_indels_insertion
