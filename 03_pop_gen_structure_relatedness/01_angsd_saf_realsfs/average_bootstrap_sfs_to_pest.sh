###
# average_bootstrap_sfs_to_pest.sh does this:
# Averages multiple SFS rows into a single mean prior vector suitable for use as a pest file.
###
awk '{
  for(i=1;i<=NF;i++) s[i]+=$i
  n++
}
END{
  for(i=1;i<=NF;i++){
    printf "%.10f%s", s[i]/n, (i<NF?OFS:ORS)
  }
}' OFS=" " catfish.global.folded.sfs > catfish.global.folded.sfs.mean.pest
