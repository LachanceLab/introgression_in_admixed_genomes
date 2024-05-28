#!/usr/bin/bash
aa_only=false
header=false
while getopts "b:i:g:e:c:o:s:ah" arg; do
  case $arg in
    b)
      bedtools_path=$OPTARG;;
    i)
      ibdmix=$OPTARG;;
    g)
      genes=$OPTARG;;
    e)
      ensembl=$OPTARG;;
    c)
      genomefile=$OPTARG;;
    o)
      output=$OPTARG;;
    s)
      slop=$OPTARG;;
    a)
      aa_only=true;;
    h)
      header=true;;
  esac
done
if [[ ${header} == true ]]; then
  if [[ ${aa_only} == true ]]; then
    ${bedtools_path} intersect -a <(${bedtools_path} merge -i <(grep -w AA ${ibdmix} | sort -k1,1 -k2,2n | cut -f1-3 | ${bedtools_path} slop -b ${slop} -i - -g ${genomefile})) -b ${genes} -wao | awk -F '\t' '{if ($7 != ".") print $7}' | sort | uniq | xargs -I {} grep -w {} ${ensembl} | cut -f2 | sort | uniq > ${output}
  else
    ${bedtools_path} intersect -a <(${bedtools_path} merge -i <(sort -k1,1 -k2,2n ${ibdmix} | cut -f1-3 | head -n -1 | ${bedtools_path} slop -b ${slop} -i - -g ${genomefile})) -b ${genes} -wao | awk -F '\t' '{if ($7 != ".") print $7}' | sort | uniq | xargs -I {} grep -w {} ${ensembl} | cut -f2 | sort | uniq > ${output}
  fi
else
  if [[ ${aa_only} == true ]]; then
    ${bedtools_path} intersect -a <(${bedtools_path} merge -i <(grep -w AA ${ibdmix} | sort -k1,1 -k2,2n | cut -f1-3 | ${bedtools_path} slop -b ${slop} -i - -g ${genomefile})) -b ${genes} -wao | awk -F '\t' '{if ($7 != ".") print $7}' | sort | uniq | xargs -I {} grep -w {} ${ensembl} | cut -f2 | sort | uniq > ${output}
  else
    ${bedtools_path} intersect -a <(${bedtools_path} merge -i <(sort -k1,1 -k2,2n ${ibdmix} | cut -f1-3 | ${bedtools_path} slop -b ${slop} -i - -g ${genomefile})) -b ${genes} -wao | awk -F '\t' '{if ($7 != ".") print $7}' | sort | uniq | xargs -I {} grep -w {} ${ensembl} | cut -f2 | sort | uniq > ${output}
  fi
fi