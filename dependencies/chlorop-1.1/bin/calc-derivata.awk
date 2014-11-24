# FILE calc-derivata.awk
#
# Calculates a numerical derivative of the S-score for how-test.out -files.
#
# d is the size of one HALF of the window within which the derivative
# is calculated. 
# (The position "i" corresponds to the assignment "C" in how-files.)
# AUTHOR: Olof Emanuelsson, olof@sbc.su.se

BEGIN {
  out_test=0
  if (d==0) {
      exit
      print "exit"
  }
}

{ out_single=0 }

/^ # / {
  name =$2
  len = $3
  Nseq++
  max_len=83           # search highest derivative within # aa; 83 => CS located within 100 first aa
  highest=0
  highest_aa=0
  if (len<max_len) max_len=len
  getline

# Read the NN (HOW) scores
  for (i=1; i<= len; i++) {
          aa_nr[i]=$1 
	  aa[i]=$2
	  nn_score[i]=$5
	  if (i<len) getline	  
  }

# calculate the derivative and get the position and value of the highest derivative
  for (i=1; i <= len; i++) {
    first_sum=0
    second_sum=0
    if (i>d && i<=len-d+1) {
      for (j=i-d; j<=i-1; j++) {
	first_sum += nn_score[j]
      }

      for (j=i; j<=i+d-1; j++) {
	second_sum += nn_score[j]
      }
    }
    if (i<=d || i>len-d+1){
      derivata[i]= 0
    }
    else {
      derivata[i]= (first_sum-second_sum)/d
    }
  }

# search highest derivative (within the limit set by "max_len")
  for (i=1; i <= max_len; i++) {
    if (derivata[i] > highest) { 
      highest=derivata[i]
      highest_aa= i
    }     
  }

# print results 
  printf "\n# %s %5d %5.3f %5d \n", name, len, highest, highest_aa
  for (i=1; i<= len; i++) {
    printf "%4d %1s %5.3f %5.3f\n", aa_nr[i], aa[i], nn_score[i], derivata[i]
  }
}

