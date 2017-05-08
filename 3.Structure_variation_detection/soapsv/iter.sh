# argv 1 input file
# argv 2 fa_stat file
# argv 3 fix name
# argv 4 query seq file

i=95
#echo "axtdetail2axt_nt $1 95 >all_${3}_i$i"
echo "./axt_correction $1 > ${1}_corrected"
echo "./lastz_pipeline/axtSort ${1}_corrected all_${3}_i$i""_sort"
echo "./lastz_pipeline/axtBest all_${3}_i$i""_sort all all_${3}_i$i""_sort_best"
echo "./lastz_postprocess all_${3}_i$i""_sort_best 3 1>all_${3}_i$i""_sort_best_pp 2>/dev/null"
echo "./trash_comments all_${3}_i$i""_sort_best_pp >all_${3}_i$i""_sort_best_pp_nc"
echo "./trans_pos.pl $2 all_${3}_i$i""_sort_best_pp_nc >all_${3}_i$i""_sort_best_pp_nc_tp"
echo "./msort -L4 -km5,n6,rn7 all_${3}_i$i""_sort_best_pp_nc_tp >all_${3}_i$i""_sort_best_pp_nc_tp_ms"
echo "./best_hit all_${3}_i$i""_sort_best_pp_nc_tp_ms >all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh"
echo "./best_hit_2 all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh >all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_bh2"
echo "awk 'NF==9' all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh >all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_oi"
echo "awk 'NF==9' all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_bh2 >all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_bh2_oi"
echo "./intro_indel ${4} all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_bh2 > ${3}_intro_indels"
echo "./detect_inv.pl all_${3}_i$i""_sort_best_pp_nc_tp_ms_bh_bh2_oi ${3}_inversion"
