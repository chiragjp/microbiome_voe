#config file
IBD	FALSE	TRUE
T2D	FALSE	TRUE
T2D	TRUE	TRUE
ACVD	FALSE	TRUE
cirrhosis	FALSE	TRUE
adenoma	TRUE	TRUE
adenoma	FALSE	TRUE
CRC	TRUE	TRUE
CRC	FALSE	TRUE
otitis	FALSE	TRUE
IBD	FALSE	FALSE
T2D	FALSE	FALSE
T2D	TRUE	FALSE
ACVD	FALSE	FALSE
cirrhosis	FALSE	FALSE
adenoma	TRUE	FALSE
adenoma	FALSE	FALSE
CRC	TRUE	FALSE
CRC	FALSE	FALSE
otitis	FALSE	FALSE




#submit initial regressions

for file in t_*.rds; do while read p; do echo sbatch -p short -c 1 -n 1 --mem=10G -t 0-00:15 ./runr.sh univariate_regression.R $file $p; done<phenotypes_randomize_batch; done

#build config files for ml
while read p; do

while read pp; do

while read ppp; do

ls modeling*"${p}"*"${pp}"*"${ppp}"*rds | grep -v overlapping > ml_config_file_"${p}"_"${pp}"_"${ppp}"

done<baz

done<bar

done<foo

while read p; do

while read pp; do


ls modeling*"${p}"*"${pp}"*overlapping*rds | grep overlapping > ml_config_file_"${p}"_"${pp}"_overlapping_genefamilies.rds


done<bar

done<foo

for file in ml_config_file*metaphlan;

do

pheno=$(echo $file | cut -f4 -d_)

echo sbatch -p short -c 1 -n  1 --mem=15G -t 0-07:00 ./runr.sh ml_analysis.R $file t_metaphlan_full.rds $pheno >> fulljoblist

done

for file in ml_config_file*pathway;

do

pheno=$(echo $file | cut -f4 -d_)

echo sbatch -p short -c 1 -n  1 --mem=120G -t 0-11:00 ./runr.sh ml_analysis.R $file t_pathways_full.rds $pheno >> fulljoblist

done

for file in $(ls ml_config_file*genefam* | grep -v overlapping);

do

pheno=$(echo $file | cut -f4 -d_)

echo sbatch -p short -c 1 -n  1 --mem=230G -t 0-11:00 ./runr.sh ml_analysis.R $file t_genefamilies_full.rds $pheno >> fulljoblist

done
