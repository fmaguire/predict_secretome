# FILE analys2data_A.awk
#
# IN:	file with averaged HOW network output (results/5to1HOW_scores.$$)
# OUT:	file with HOWLIN-input entries (only numbers...) AND
#	file with name, length and no of sequences
# NOTE:	if an input sequence is less than 100 aa long, the lacking stretch
# will be filled in with 0.000 scores. In the training and testing of the 
# howlin network in the postprocessing, this is the case only with one 
# sequence (cTP).
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se

BEGIN {
	counter=0
	max_name_length=20
}

/^ #/ {
	counter++
	name = $2
	len = $3
	name_toprint=name
	for (i=length(name); i<max_name_length; i++) {
		name_toprint = name_toprint " "
	}
	printf "%20s   %5.0f\n", name_toprint, len >> name_file
	i=0
	cTP_transit=0
	no_transit=0
	while  (i<len) {
		getline
		i++
		if (i<=w) printf "%5.3f  ", $5			
	}
	if (len<w) {
		for (j=1; j<=w-len; j++) {		  # if <100, zeroes filled in
			if (i+j<=w) printf "0.000  "
		}
	}	
	printf "%2.0f %2.0f\n",cTP_transit,no_transit
	getline		# to get it right with the 5to1HOW_scores.$$-file
	next
}
