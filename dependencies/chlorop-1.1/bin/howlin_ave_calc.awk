# FILE howlin_ave_calc.awk
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se

BEGIN {		
		cut_off=0.5
                for (i=1; i<=noofseq; i++) score[i]=0
}
        /^ #/ { 
                for (i=0; i<=4 ; i++) {
                        score[$2]+=$(i*8+6)
                }
		# '/=' does not work with certain AWK's
                score[$2] = score[$2]/5
        }

END {
	for (i=1; i<=noofseq; i++) { 
		if (score[i] > cut_off) printf "   %5.3f   Y\n", score[i]
		else 		    	printf "   %5.3f   -\n", score[i]
	}							
}
