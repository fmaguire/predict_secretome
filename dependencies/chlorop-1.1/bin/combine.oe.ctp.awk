# FILE combine.oe.ctp.awk
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se

/^ #/ { 
	N++
	len = substr($0,39,6)+0 
	name = substr($0,7,20)
	seq = ""
}

/SINGLE/ { 
	out=1 
	next
}

(out) {
	S[$1] = 0
	AA[$1]= $2
	if ($3== "P" || $3=="_")  TP_real[$1]=$3
	else TP_real[$1]="?"
	for (i=0; i<5; i++) {
		S[$1] += $(6*i+5)
	}
	# '/=' does not work with certain AWK's
	S[$1] = S[$1]/5
	if (S[$1]>1) S[$1]=1
}

out && $1 == len { 
	out=0
	printf " # %-10s %5d \n",
		name, len 
	for (i=1; i<=len; i++){
	        if (S[i]>= 0.5) TP_assign = "P"
		else TP_assign = "_"
	        printf "%5d %-1s %-1s %-1s %5.3f \n", i, AA[i], TP_real[i], TP_assign, S[i] 
		# printing of TP_real lacks of course meaning when using unknown sequences
	}
	printf "\n" 
}
