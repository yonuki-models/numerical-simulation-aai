#!/bin/sh

var='T'
base_dir='../../anime/exRo10N10/'
tmp_dir=${base_dir}'tmp/'
mkdir -p ${tmp_dir}

for i in `seq 0 100`
do
  ip=$(($i + 200))
  si=$(printf "%04d\n" "${i}")
  sip=$(printf "%04d\n" "${ip}")
  cp ${base_dir}${var}${sip}.png ${tmp_dir}${var}${si}.png
  echo ${i}
done
ffmpeg -framerate 10 -i ${tmp_dir}${var}%04d.png -vcodec libx264 -pix_fmt yuv420p -r 60 ${var}.mp4

rm ${tmp_dir}${var}*
