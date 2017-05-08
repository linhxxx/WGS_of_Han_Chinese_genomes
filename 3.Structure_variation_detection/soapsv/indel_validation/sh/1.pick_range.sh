#argv1 extending length
#argv2 SV list
#argv3 inversion list
echo "awk '"'$2=="Insertion"{start=$3-'"$1"';if(start<1){start=1};print $1","$3","$4"\t"$1","start","$4+'"$1"';}'"' $2 >insertion_scaf.origin"
echo "awk '"'$2=="Deletion"{start=$3-'"$1"';if(start<1){start=1};print $1","$3","$4"\t"$1","start","$4+'"$1"';}'"' $2 >deletion_scaf.origin"
echo "awk '"'$2=="Insertion"{start=$6-'"$1"';if(start<1){start=1};print $5","$6","$7"\t"$5","start","$7+'"$1"';}'"' $2 >insertion_hg18.origin"
echo "awk '"'$2=="Deletion"{start=$6-'"$1"';if(start<1){start=1};print $5","$6","$7"\t"$5","start","$7+'"$1"';}'"' $2 >deletion_hg18.origin"
echo "awk '"'{start=$3-'"$1"';if(start<1){start=1};print $1","$3","$4"\t"$1","start","$4+'"$1"';}'"' $3 >inversion_scaf.origin"
echo "awk '"'{start=$6-'"$1"';if(start<1){start=1};print $5","$6","$7"\t"$5","start","$7+'"$1"';}'"' $3 >inversion_hg18.origin"
echo "awk '"'{print $2}'"' insertion_scaf.origin >insertion_scaf.list"
echo "awk '"'{print $2}'"' deletion_scaf.origin >deletion_scaf.list"
echo "awk '"'{print $2}'"' insertion_hg18.origin >insertion_hg18.list"
echo "awk '"'{print $2}'"' deletion_hg18.origin >deletion_hg18.list"
echo "awk '"'{print $2}'"' inversion_scaf.origin >inversion_scaf.list"
echo "awk '"'{print $2}'"' inversion_hg18.origin >inversion_hg18.list"
