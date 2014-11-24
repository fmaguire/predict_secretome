# FILE calc_of_motif_score.awk

# Reads matrix file, short deriv-list file and calculates scores for every 
# sequence and then determines the best score and presents the corresponding
# sequence and cleavage site.
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se

BEGIN {
	Nseq=0
}

/^##/ {
    MEME_alength=20 
    win_len=$3
    cs_lok=$5 		# The cs position in window with the motif found by MEME
    for (i=1; i<=win_len; i++) {
        getline
	for (j=1; j<=MEME_alength; j++ ) {
	  MEME_matrix[(i-1)*MEME_alength + j]= $(j+1)
	}
    }
}


/^# / {	# Calculate the scores and search the highest
	Nseq++
	highest_score[Nseq]=-10000 # These 2 variables are used in the search 
	highest_start_pos[Nseq]=0  # for the highest SCORE (from the sc. matrix)
	name[Nseq]=$2
	highest_deriv_aa[Nseq]=$3
	getline
	sequence[Nseq]= $1
	seq_len[Nseq]=length(sequence[Nseq])
	
	for (i=1; i<=seq_len[Nseq]; i++) {
		score[i]=0
	}
	
		# Prepare file with cleavage site scores
	print "\n# " name[Nseq] >> csscorefile  # will contain rediduewise cs-score
	printf "   1   --- \n   2   --- \n" >> csscorefile
		# Calculate scores
	for (i=1; i<=seq_len[Nseq]; i++) {
		for (j=1; j<=win_len; j++) {
			score[i]=score[i]+get_score( substr(sequence[Nseq],i+j-1,1)  ,j)
		}
		if (i<=seq_len[Nseq]-2) printf "%4d %6.3f\n",i+2,score[i] >> csscorefile
	}

		# Search highest score in area +-20 aa around highest derivative
	if (highest_deriv_aa[Nseq]>win_in)  start_offset=win_in	# getting the starting aa
	else start_offset = highest_deriv_aa[Nseq]-1	
	for (i=highest_deriv_aa[Nseq]-start_offset; i<=highest_deriv_aa[Nseq]+win_in-1; i++) {
		if(score[i]> highest_score[Nseq]){
			highest_score[Nseq]=score[i]
			highest_start_pos[Nseq]=i
		}
	}
	
		# Calculate the cTP len predicted by this method
	tp_len_assign[Nseq]=highest_start_pos[Nseq] + cs_lok - 1
}


function get_score(AA, position,	aa_type) {
	if (AA=="A") aa_type=1
	if (AA=="C") aa_type=2
	if (AA=="D") aa_type=3
	if (AA=="E") aa_type=4
	if (AA=="F") aa_type=5
	if (AA=="G") aa_type=6
	if (AA=="H") aa_type=7
	if (AA=="I") aa_type=8
	if (AA=="K") aa_type=9
	if (AA=="L") aa_type=10
	if (AA=="M") aa_type=11
	if (AA=="N") aa_type=12
	if (AA=="P") aa_type=13
	if (AA=="Q") aa_type=14
	if (AA=="R") aa_type=15
	if (AA=="S") aa_type=16
	if (AA=="T") aa_type=17
	if (AA=="V") aa_type=18
	if (AA=="W") aa_type=19
	if (AA=="Y") aa_type=20
	if (AA=="-") return -100 # Penalty if trying to include non-existing aa
	return MEME_matrix[(position-1)*MEME_alength + aa_type]
}

END {
	print "Nseq: " Nseq
	for (i=1; i<=Nseq; i++) {
		printf "# %10s %7.3f  Assign. cTP-len: %3d\n",name[i],highest_score[i],\
			tp_len_assign[i]
		print substr(sequence[i], highest_start_pos[i], cs_lok) " " \
			substr(sequence[i],highest_start_pos[i]+cs_lok, win_len-cs_lok) 
		print " "
	}
}
