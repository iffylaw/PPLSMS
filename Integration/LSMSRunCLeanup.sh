#!/bin/bash
# name of experiment
exp_id=$(head -1 general_settings.txt | tail -1 | awk '{print $1}')
# forcing version used (TIMESTEP or DAILY)  
forcing_frequency=$(head -1 general_settings.txt | tail -1 | awk '{print $5}')
# number of sites in experiment
NSITE=$(wc -l site_list.txt | cut -c1-8 | awk '{print $1}')

case "$forcing_frequency" in
    TIMESTEP)
	ftype=halfhourly
    ;;
    DAILY)
        ftype=daily
    ;;
    *)
    echo invalid forcing frequency selected, only TIMESTEP and DAILY allowed
    return
    ;;
esac

if [ ! -d stats ] ; then
      mkdir -p stats
fi  

outeo=${exp_id}_${ftype}_out
if [ ! -d $outeo ] ; then
      mkdir -p $outeo
fi 
mv run_jsb_site.*.? $outeo

VAR_list="PPFD Rn Qh Qle gscan GPP Reco NEE fapar"
FREQ_list="halfhourly daily monthly"
for var in $VAR_list ; do
  echo $var

  for freq in $FREQ_list ; do
    outfile=stats/${var}.${exp_id}.${ftype}.stats.${freq}.txt

    i=2
    while [ $i -le $NSITE ] ; do
      SITE=$(head -${i} site_list.txt | tail -1 | awk '{print $1}')
      if [ $i -eq 2 ] ; then
        h1=$(grep RMSE ${SITE}/Output/${SITE}.${exp_id}.${ftype}.stats.${freq}.txt)
        h2=$(grep JSBACH ${SITE}/Output/${SITE}.${exp_id}.${ftype}.summary.txt)
        echo \"SITE\" $h1 $h2 > ${outfile}
      fi 
      d=$(grep ${var} ${SITE}/Output/${SITE}.${exp_id}.${ftype}.stats.${freq}.txt)
      e=$(grep ${var} ${SITE}/Output/${SITE}.${exp_id}.${ftype}.summary.txt | cut -d ' ' -f2-4)
      echo \"${SITE}\" $d $e >> ${outfile}
      let i=i+1
    done
  done
done
