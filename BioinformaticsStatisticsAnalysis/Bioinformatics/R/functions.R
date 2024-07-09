aa_list = c('A', 'C', 'D', 'E', 'F','G', 'H','I', 'K', 'L', 'M','N','P', 'Q', 'R', 'S', 'T','V','W', 'Y', '*')
aa_freq = c( 4,   2,   2,   2,   2,  4,   2,  3,   2,   6,   1,  2,  4,   2,   6,   6,   4,  4,  1,   2,   3)
aa_freq =  aa_freq/64

convert_3_to_1 = function(vec){
  new_vec = gsub("Pro", "P", vec)
  new_vec = gsub("Ser", "S", new_vec)
  new_vec = gsub("Phe", "F", new_vec)
  new_vec = gsub("Ala", "A", new_vec)
  new_vec = gsub("Ile", "I", new_vec)
  new_vec = gsub("Gln", "Q", new_vec)
  new_vec = gsub("Asp", "D", new_vec)
  new_vec = gsub("Asn", "N", new_vec)
  new_vec = gsub("Tyr", "Y", new_vec)
  new_vec = gsub("Trp", "W", new_vec)
  new_vec = gsub("Gly", "G", new_vec)
  new_vec = gsub("Met", "M", new_vec)
  new_vec = gsub("Leu", "L", new_vec)
  new_vec = gsub("His", "H", new_vec)
  new_vec = gsub("Arg", "R", new_vec)
  new_vec = gsub("Glu", "E", new_vec)
  new_vec = gsub("Val", "V", new_vec)
  new_vec = gsub("Thr", "T", new_vec)
  new_vec = gsub("Cys", "C", new_vec)
  new_vec = gsub("Lys", "K", new_vec)
  return(new_vec)
}

