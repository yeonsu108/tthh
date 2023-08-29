targetpath="./samples1/"
mkdir -p $targetpath

root -l -b -q ana_tthh_df.C'("Semi_tthh", "'$targetpath'")' &> ./$targetpath/log_tthh &
root -l -b -q ana_tthh_df.C'("Semi_ttbb", "'$targetpath'")' &> ./$targetpath/log_ttbb &
root -l -b -q ana_tthh_df.C'("Semi_ttbbbb", "'$targetpath'")' &> ./$targetpath/log_ttbbbb &
root -l -b -q ana_tthh_df.C'("Semi_ttbbcc", "'$targetpath'")' &> ./$targetpath/log_ttbbcc &

#targetpath="./"
#mkdir -p $targetpath
#root -l -b -q ana_tthh_df.C'("Test_tthh", "'$targetpath'")' &> ./$targetpath/log_tthh &
