# argv1 input scaffold fa
# argv2 selected list directory
# argv3 output directory
# argv4 task name

echo "cp /panfs/RD/luoruibang/bin/lastz_pipeline/2.pick_novel/* ."
echo "cd $3"
echo 'grep ">" '"$1"" | perl -ne '"'chomp; />(\S+)/; print "$1\n";'"' >all.list";
echo "cat $2/chr* >selected.list"
echo "sort selected.list | uniq >selected.uniq.list"
echo "./pick_novelscaf all.list selected.uniq.list >novel.list"
echo "/panfs/RD/luoruibang/bin/lastz_pipeline/1.blat/list/pick_scaf novel.list $1 novel.fa"
echo "/panfs/RD/luoruibang/bin/fasta_seperator novel.fa nvs ./ 5000000"
echo "ls $3/nvs* > fa.list"
echo "mkdir result"
echo "./scriptmaker /panfs/RD/luoruibang/temp/reference/hg18/capsule/capsule19.list fa.list $3/result >run.sh"
echo "runcount=\`wc -l run.sh\`"
#echo "job -i run.sh -c "'$runcount'" -n $4 -r 1.9"
