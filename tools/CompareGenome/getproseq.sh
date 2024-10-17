full=$1
short=${full:0:1}.${full#*_}
abb=${short:0:1}${short:2:2}

seqkit grep -f $short/$short.final.list db_dir/$full/$short.pep > $short/${abb}GIS3.final.pep
