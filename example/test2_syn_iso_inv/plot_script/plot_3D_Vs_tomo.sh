#!/bin/bash
# script for ploting 3D Vs slices
# Author: Chuanming Liu (USTC)

cpt_path=.
input_path=.
output_path=.
tomo_cpt_in=$cpt_path/"BlueDarkRed18.cpt"
tomo_cpt='blue_red.cpt'

J=M6.5c
R_tomo=101.25/104.75/23/26.5

azm_file[1]=$input_path/DSurfTomo.inv
azm_file[2]=$input_path/DSurfTomo.true
fileNum=2
# depth
depth[1]=0
depth[2]=10
depth[3]=35
depth[4]=60
Tnum=3

offset[1]="-Y5c -X2.5c"
offset[2]=-X8.2c
offset[3]=-X8.2c
offset[4]=-X8.2c
offset[5]="-X-24.6c -Y-11.3c"
offset[6]="-X8.2c"
offset[7]="-X8.2c"
offset[8]="-X8.2c"

# psscale
R_T=2.8/4.2/0.02
psscale_title="Vs (km/s)"
psscale_B=xa0.4f0.1
psscale_By=y+l"$psscale_title"

# GMT 5 defaults
gmt defaults -D > .gmtdefaults4
gmt set PS_MEDIA A4
gmt set MAP_FRAME_TYPE Plain
gmt set MAP_FRAME_PEN 2p
gmt set MAP_FRAME_WIDTH 0.11c
gmt set FONT Times-Roman
gmt set FONT_ANNOT_PRIMARY 21p,Times-Roman 
gmt set FONT_TITLE 22p,Times-Roman
gmt set FONT_LABEL 19p,Times-Roman
gmt set MAP_TITLE_OFFSET 0.0c
gmt set MAP_GRID_PEN_PRIMARY thinner,black
gmt set COLOR_BACKGROUND 0/0/255
gmt set COLOR_FOREGROUND 255/0/0
gmt set COLOR_NAN 255/255/255

for ((ff=1; ff<=$fileNum; ff=ff+1)); do
azmthfile=${azm_file[$ff]}
pic_name=$output_path/${azm_file[$ff]}.ps
echo ${azm_file[$ff]}

gmt makecpt -C$tomo_cpt_in -T$R_T -I  > $tomo_cpt
gmt psxy -J$J -R$R_tomo -T  -K > $pic_name
for ((i=1; i<=$Tnum; i=i+1)); do
    title="Depth "${depth[$i]}" km"
    echo $title
    awk '{if($3==depth1) print $1,$2,$4}' depth1=${depth[$i]} $azmthfile| gmt surface -R$R_tomo -I0.01  -Gtomo_grd
    gmt grdimage tomo_grd -J$J -R$R_tomo -C$tomo_cpt  -Bxa1f1 -Bya1f1 -BWeSn+t"$title" ${offset[$i]} -K  -O >> $pic_name
    if  [ $i = 2 ]  ;then
    gmt psscale -Dx-0.5c/-2c+w7c/0.4c+h+e -C$tomo_cpt -B$psscale_B -B"$psscale_By" -K -O >> $pic_name
    fi
    rm tomo_grd
done
gmt psxy  -J$J -R$R_tomo -T -O  >> $pic_name
gmt psconvert -P -A1c -Tf  $pic_name
done
echo "finish plot."

rm $input_path/gmt.history $input_path/gmt.conf
rm $input_path/*.ps
rm $input_path/$tomo_cpt
