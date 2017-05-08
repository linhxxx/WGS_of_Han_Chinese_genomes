#argv1 s_avg_ref
#argv2 p_avg_ref
#argv3 s_avg_scaf
#argv4 p_avg_scaf
#argv5 depth_limit

scrPATH/indel_validation/indel_validation.pl insertion insertion.combined insertion.validated scaf 50 $1 $2 $3 $4 $5
scrPATH/indel_validation/indel_validation.pl deletion deletion.combined deletion.validated ref 50 $1 $2 $3 $4 $5
scrPATH/indel_validation/indel_validation.pl inversion inversion.combined inversion.validated scaf 50 $1 $2 $3 $4 $5
