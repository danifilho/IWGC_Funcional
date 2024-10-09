#!/bin/bash
input_list=$1
fasta=$2
go_list=$3
xdir=$4
volume_name=$5

for acc in `cat ${input_list}  `
   do
       #this should put work into the background, but we have to wait for it all
       # to finish so
       accession=$( printf "${acc}" | sed 's/-/./g' )
       #adjust the directories
       singularity exec -B ${volume_name}:/data $volume_name/images/samtools.sif /bin/bash -c "cd /data && samtools faidx /data/${fasta} \"${acc}\" > /data/splits/${xdir}/${accession}.fa"
       printf "${acc}\n" >${accession}.txt 
       #this should put work into the background, but we have to wait for it all
       # to finish so
       egrep -wf  ${accession}.txt $volume_name/${go_list} > ${xdir}/${accession}.goterms
       rm ${accession}.txt

    if [ ! -d "${xdir}/tmp_${accession}" ] ; then mkdir "${xdir}/tmp_${accession}" ; fi

        singularity exec \
        --bind ${volume_name}/splits:/input,${xdir}/tmp_${accession}:/tmp \
        ${volume_name}/images/multiloc2_v3.img \
        python2 /opt/MultiLoc2-1/MultiLoc2/src/multiloc2_prediction.py \
              -fasta=/input/${xdir}/${accession}.fa \
              -output=advanced \
              -go=/input/${xdir}/${accession}.goterms \
              -origin=plant \
              -result=/input/${accession}.multLoc2out >> ${xdir}/log.txt 2>> ${xdir}/error.txt      

    if [ -f ${accession}.multLoc2out  ] ; then
      cat ${accession}.multLoc2out | awk ' ! a[$1]++ ' \
      | egrep -v '^Detailed|^sequence|^predictor|^origin|^$|^Multi' >> ${xdir}/ml2.out
    else
      echo ${accession} >> ${xdir}/failed.list
    fi
   rm -rf ${xdir}/${accession}.fa ${xdir}/${accession}.goterms ${xdir}/tmp_${accession}
 done
