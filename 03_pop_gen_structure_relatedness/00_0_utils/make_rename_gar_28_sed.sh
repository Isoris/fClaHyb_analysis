###
# make_rename_gar_28_sed.sh does this:
# Generates a sed substitution file to rename C. gariepinus chromosome headers
# from accession-based names to short LG names (LG01-LG28).
###
cat > rename_gar_28.sed <<'SED'
s/^C_gar_CM107359\.1_LG01\b/C_gar_LG01/
s/^C_gar_CM107360\.1_LG02\b/C_gar_LG02/
s/^C_gar_CM107361\.1_LG03\b/C_gar_LG03/
s/^C_gar_CM107362\.1_LG04\b/C_gar_LG04/
s/^C_gar_CM107363\.1_LG05\b/C_gar_LG05/
s/^C_gar_CM107364\.1_LG06\b/C_gar_LG06/
s/^C_gar_CM107365\.1_LG07\b/C_gar_LG07/
s/^C_gar_CM107366\.1_LG08\b/C_gar_LG08/
s/^C_gar_CM107367\.1_LG09\b/C_gar_LG09/
s/^C_gar_CM107368\.1_LG10\b/C_gar_LG10/
s/^C_gar_CM107369\.1_LG11\b/C_gar_LG11/
s/^C_gar_CM107370\.1_LG12\b/C_gar_LG12/
s/^C_gar_CM107371\.1_LG13\b/C_gar_LG13/
s/^C_gar_CM107372\.1_LG14\b/C_gar_LG14/
s/^C_gar_CM107373\.1_LG15\b/C_gar_LG15/
s/^C_gar_CM107374\.1_LG16\b/C_gar_LG16/
s/^C_gar_CM107375\.1_LG17\b/C_gar_LG17/
s/^C_gar_CM107376\.1_LG18\b/C_gar_LG18/
s/^C_gar_CM107377\.1_LG19\b/C_gar_LG19/
s/^C_gar_CM107378\.1_LG20\b/C_gar_LG20/
s/^C_gar_CM107379\.1_LG21\b/C_gar_LG21/
s/^C_gar_CM107380\.1_LG22\b/C_gar_LG22/
s/^C_gar_CM107381\.1_LG23\b/C_gar_LG23/
s/^C_gar_CM107382\.1_LG24\b/C_gar_LG24/
s/^C_gar_CM107383\.1_LG25\b/C_gar_LG25/
s/^C_gar_CM107384\.1_LG26\b/C_gar_LG26/
s/^C_gar_CM107385\.1_LG27\b/C_gar_LG27/
s/^C_gar_CM107386\.1_LG28\b/C_gar_LG28/
SED
