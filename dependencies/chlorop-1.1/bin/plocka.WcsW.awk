# FILE plocka.WcsW.awk
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se


/^ [0-9]* / {
    len=$1
    name=substr($2,1,11)
    deriv_max=predicted_tp_len(name) 
    getline
    sekv=$1
    temp_len=len
    while (temp_len/80 > 1) {
	getline
	sekv= sekv $1
	temp_len=temp_len-80
    }
    # Print out HOW-style output
    printf "# %s %3d\n",name,deriv_max
    print sekv
}

function predicted_tp_len(name) {
    getline<deriv_file
    deriv_max=$5
    return deriv_max
}

