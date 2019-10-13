#!/bin/bash
mkdir -p replicated_dvg
mkdir -p replicated_dvg/ind
mkdir -p replicated_dvg/final
mkdir -p replicated_dvg/depth
while getopts ":1:2:" opt; do
  case "$opt" in
    1 ) one=$OPTARG      ;;
    2 ) two=$OPTARG     ;;
esac
done
shift $(( OPTIND - 1 ))
replicates=($one $two)
# Removes old files if they exist
rm -f replicated_dvg/total_average_depth.tsv
rm -f replicated_dvg/all_replicated.tsv

for i in ${one}/dvg_summary/raw/*.tsv; do
samples+=($(echo $i | rev | cut -d '/' -f 1 | rev | cut -d'.' -f1))
done

for i in ${replicates[@]}; do
out_name=$( echo $i | sed 's/\//_/g' )
# Combines output files from all replicates
cat $i/dvg_summary/raw/*.tsv > replicated_dvg/all_${i}_dvg.tsv

awk '{OFS="\t"; print $1,$4,$5,$6}' replicated_dvg/all_${i}_dvg.tsv | \
sort > replicated_dvg/all_${i}_dvg_clean.tsv
done

# Finds DVGs (sample_day_dvg) which are present in both replicates
comm -12 replicated_dvg/all_${replicates[0]}_dvg_clean.tsv \
replicated_dvg/all_${replicates[1]}_dvg_clean.tsv | \
sed -e 's/_/\t/g' - > replicated_dvg/replicated_dvgs.tsv

for i in ${samples[@]}; do
    counter=0
    cmd="paste "
    # Greps list of replicated DVGs against the original DVG lists
    for j in ${replicates[@]}; do
        grep -f replicated_dvg/replicated_dvgs.tsv ${j}/dvg_summary/raw/${i}.norm.dvgCount.tsv > \
        replicated_dvg/ind/${i}.${j}.replicated.tsv
    done
    # Combines the two lists of repeated DVGs, summing read support and averaging relative read support
    for j in ${replicates[@]}; do 
        if [[ "$counter" = 0 ]]; then
            cmd+="<(sort replicated_dvg/ind/${i}.${j}.replicated.tsv) "
            counter=$((counter+1))
        else
            cmd+="<(sort replicated_dvg/ind/${i}.${j}.replicated.tsv | awk '{print \$7,\$8}') "
        fi
    done
    cmd+=" | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7+\$9,(\$8+\$10)/2}' > \
    replicated_dvg/final/${i}.final.replicated.tsv"
    eval "$cmd"
done

# Sums by position read depth from individual replicates
for i in ${samples[@]}; do
    counter=0
    cmd="paste "
    for j in ${replicates[@]}; do
        if [[ "$counter" = 0 ]]; then
            cmd+="<(cat ${j}/depth/${i}.all.depth.tsv) "
            counter=$((counter+1))
        else
            cmd+="<(awk '{print \$3,\$4}' ${j}/depth/${i}.all.depth.tsv) "
        fi
    done
    cmd+=" | awk '{OFS=\"\t\"; print \$1,\$2,\$3+\$5,\$4+\$6}' > \
    replicated_dvg/depth/${i}.depth.tsv"
    eval "$cmd"
done

# Tabulates read depth statistics
for i in ${samples[@]}; do 
echo $i
awk -v subject_day="${i}" \
'{g[1]="WG";
g[2]="PB2"; 
g[3]="PB1";
g[4]="PA";
g[5]="HA";
g[6]="NP";
g[7]="NA";
g[8]="M";
g[9]="NS"}{split($1,a,"|"); \
b["WG"]+=1; c["WG"]+=$3; d["WG"]+=$3*3; e["WG"]+=$3; f["WG"]+=$4*$4; \
b[a[1]]+1; c[a[1]]+=$3; d[a[1]]+=$3*$3; e[a[1]]+=$4; f[a[1]]+=$4*$4}END\
{OFS=""; for (i in g){print subject_day}}' 
done

for i in ${samples[@]}; do
awk -v subject_day="${i}" \
'{g[1]="WG";
g[2]="PB2"; 
g[3]="PB1";
g[4]="PA";
g[5]="HA";
g[6]="NP";
g[7]="NA";
g[8]="M";
g[9]="NS"}{split($1,a,"|"); \
b["WG"]+=1; c["WG"]+=$3; d["WG"]+=$3*$3; e["WG"]+=$4; f["WG"]+=$4*$4; \
b[a[1]]+=1; c[a[1]]+=$3; d[a[1]]+=$3*$3; e[a[1]]+=$4; f[a[1]]+=$4*$4}END\
{OFS=""; for (i in g){print \
subject_day,"\t",g[i], "\t",c[g[i]]/b[g[i]], "\t", sqrt (d[g[i]]/b[g[i]] - (c[g[i]]/b[g[i]])^2), "\t", 
e[g[i]]/b[g[i]], "\t", sqrt (f[g[i]]/b[g[i]] - (e[g[i]]/b[g[i]])^2), "\n"}}' \
replicated_dvg/depth/${i}.depth.tsv | sed '/^[[:space:]]*$/d' - >> \
replicated_dvg/total_average_depth.tsv
done



# Calculates average read length from the two replicates
counter=0
cmd="paste "
for i in ${replicates[@]}; do 
    if [[ "$counter" = 0 ]]; then
        cmd+="<(sort $i/read_length/all_average_read_length.tsv) "
        counter=$((counter+1))
    else
        cmd+="<(sort $i/read_length/all_average_read_length.tsv) | awk '{OFS=\"\t\"; \\
print \$1, \$2, (\$4+\$9)/(\$3+\$8), sqrt ((\$5/\$3-(\$4/\$3)^2)/\$3 + (\$10/\$8 - (\$9/\$8)^2)/\$8)}' "
    fi
done
cmd+=" > replicated_dvg/total_average_read_length.tsv"
eval "$cmd"

# Concatenates all replicated dvg data into a single file
cat replicated_dvg/final/*.tsv > replicated_dvg/all_replicated.tsv













