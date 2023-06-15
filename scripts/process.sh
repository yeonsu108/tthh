
targetpath="./samples/"
mkdir -p $targetpath

root -l -b -q ana_tthh_df.C'("tthh", "'$targetpath'")' &> ./$targetpath/log_tthh &
root -l -b -q ana_tthh_df.C'("ttbb", "'$targetpath'")' &> ./$targetpath/log_ttbb &
root -l -b -q ana_tthh_df.C'("ttbbbb", "'$targetpath'")' &> ./$targetpath/log_ttbbbb &
root -l -b -q ana_tthh_df.C'("ttbbcc", "'$targetpath'")' &> ./$targetpath/log_ttbbcc &
