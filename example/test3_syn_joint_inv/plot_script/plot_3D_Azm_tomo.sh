#!/bin/bash
# script for ploting 3D Azi slices
# Author: Chuanming Liu (USTC)

cpt_path=.
input_path=.
output_path=.
tomo_cpt=$cpt_path/"BlueDarkRed18.cpt"
cpt_in=./tomo_mk.cpt


file[1]=Gc_Gs_model.real
file[2]=Gc_Gs_model.inv
azm_file[1]=$input_path/${file[1]}
azm_file[2]=$input_path/${file[2]}
fileNum=2


# Plot setting
J=M6.5c
R_tomo=101.25/104.75/23/26.5
# depth
depth[0]=0
depth[1]=10
depth[2]=35
depth[3]=60
Tnum=3

offset[1]="-Y5c -X2c"
offset[2]=-X8.3c
offset[3]=-X8.3c
offset[4]="-X-19.6c -Y-11.3c"
offset[5]="-X8.2c"
offset[6]="-X8.2c"
offset[7]="-X8.2c"
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
gmt set FONT Times-Roman
gmt set FONT_ANNOT_PRIMARY 21p,Times-Roman 
gmt set FONT_TITLE 25p,Times-Roman
gmt set FONT_LABEL 19p,Times-Roman  
gmt set MAP_TICK_LENGTH_PRIMARY 6p/1p  
gmt set MAP_TICK_PEN_PRIMARY  2p,black
gmt set MAP_TITLE_OFFSET 0p
gmt set MAP_GRID_PEN_PRIMARY thicker,black
gmt set COLOR_BACKGROUND 0/0/255
gmt set COLOR_FOREGROUND 255/0/0
gmt set COLOR_NAN 255/255/255


for ((ff=1; ff<=$fileNum; ff=ff+1))  ;do
    azmthfile=${azm_file[$ff]}
    pic_name=${file[$ff]}.ps
    echo $azmthfile
    gmt psxy -J$J -R$R_tomo -T  -K > $pic_name
    for ((i=1; i<=$Tnum; i=i+1));do
    jj=`expr $i - 1`
    title="Depth "${depth[$jj]}'-'${depth[$i]}" km"
    echo $title

    awk '{if($3==depth1) print $1,$2,$4}' depth1=${depth[$i]} $azmthfile| gmt surface -R$R_tomo -I0.01  -Gtomo_grd

    gmt makecpt -C$tomo_cpt -T$R_T -I  >  $cpt_in

    gmt grdimage tomo_grd -J$J -R$R_tomo -C$cpt_in  -Bxa1f1 -Bya1f1 -BWeSn+t"$title" ${offset[$i]} -K  -O >> $pic_name

    awk '{if($3==depth1) print $1,$2,$5,$6*20}' depth1=${depth[$i]} $azmthfile | gmt  psxy -J$J  -R$R_tomo -SV0.2c+jc  -W1.3p,black -Gred  -K -O  >> $pic_name

    # Add lengend
    # 0.01(1%)---*40---0.1
    echo "101.69 23.2 1.5 0.7" | gmt psxy  -J$J -R$R_tomo -Sr0.4  -W0.8p,black  -Gwhite -K -O  >>$pic_name
    echo "101.45 23.2 90 0.2"  | gmt psxy  -J$J -R$R_tomo -SV0.2c+jc  -W1.3p,black -Gred  -K -O  >> $pic_name
    color=black
    echo "101.8 23.2 0 1%" | gmt pstext -R$R_tomo -J$J -F+f20p,Times--Roman,$color+a+jCM -N  -K -O  >> $pic_name

    if  [ $i = 2 ]  ;then
    gmt psscale -Dx-0.5c/-2c+w7c/0.4c+h+e -C$cpt_in -B$psscale_B -B"$psscale_By" -K -O >> $pic_name
    fi

    done
    gmt psxy  -J$J -R$R_tomo -T -O  >> $pic_name
    gmt psconvert -P -A1c -Tf  $pic_name
done

rm tomo_grd $cpt_in gmt.history gmt.conf
rm *.ps
