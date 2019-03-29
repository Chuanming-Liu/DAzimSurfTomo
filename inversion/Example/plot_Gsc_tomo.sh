#!/bin/bash
#gmt set MAP_FRAME_WIDTH 0.1c
#gmt set MAP_LABEL_OFFSET 0.1c
###############################################################
# Parameter part 1 #path
cpt_path=/Users/chuanmingliu/Lance/GMT_PIC/necessary/cpt
path=.
###############################################################
# Plot setting
J=M6c     #size for plot
tomo_cpt=$cpt_path/"BlueDarkRed18.cpt"
cpt_in=./tomo_mk.cpt
# R_tomo=101.25/104.5/23.25/26.5
R_tomo=101.4/104.4/23.49/26.4
file[1]=Gc_Gs_model.real
file[2]=Gc_Gs_model.inv
azm_file[1]=$path/${file[1]}
azm_file[2]=$path/${file[2]}
fileNum=2
# depth
depth[0]=0
depth[1]=15
depth[2]=35
depth[3]=60
Tnum=3

offset[1]="-Y5c -X1c"
offset[2]=-X8.6c
offset[3]=-X8.6c
offset[4]="-X-19.6c -Y-11.3c"
offset[5]="-X8.2c"
offset[6]="-X8.2c"
offset[7]="-X8.2c"
# psscale
R_T=2.8/4.2/0.02

psscale_title="Vs (km/s)"
psscale_B=xa0.4f0.1
psscale_By=y+l"$psscale_title"
###############################################################
# GMT 5 defaults
    gmt defaults -D > .gmtdefaults4
    gmt set PS_MEDIA A4
#   gmt set MAP_FRAME_TYPE fancy+
#   gmt set MAP_FRAME_WIDTH 2c
    gmt set MAP_FRAME_TYPE Plain
    gmt set MAP_FRAME_PEN 2p
# font
    gmt set FONT Times-Roman
    gmt set FONT_ANNOT_PRIMARY 21p,Times-Roman  # main works
    gmt set FONT_TITLE 25p,Times-Roman
    gmt set FONT_LABEL 19p,Times-Roman  # no use now
# tick
    gmt set MAP_TICK_LENGTH_PRIMARY 6p/1p  # tick length
    gmt set MAP_TICK_PEN_PRIMARY  2p,black
# offset
    gmt set MAP_TITLE_OFFSET 0p
#   gmt set MAP_GRID_PEN_PRIMARY thinner,black
    gmt set MAP_GRID_PEN_PRIMARY thicker,black
# ground color
    gmt set COLOR_BACKGROUND 0/0/255
    gmt set COLOR_FOREGROUND 255/0/0
    gmt set COLOR_NAN 255/255/255

###############################################################
# File CYCLE
   keyword="*ISO_tomo*"
   for ((ff=1;ff<=$fileNum;ff=ff+1))  ;do
       azmthfile=${azm_file[$ff]}
       pic_name=${file[$ff]}.ps
       echo $azmthfile"----------------------------"
       echo $pic_name
       gmt psxy -J$J -R$R_tomo -T -V -K > $pic_name

       for ((i=1;i<=$Tnum;i=i+1));do
             jj=`expr $i - 1`
             title="Depth "${depth[$jj]}'-'${depth[$i]}" km"
             echo $title
#---------
# GMT PLOT
             awk '{if($3==depth1) print $1,$2,$4}' depth1=${depth[$i]} $azmthfile| gmt surface -R$R_tomo -I0.01 -V -Gtomo_grd
 #           awk '{if($3==depth1) print $1,$2,$4}' depth1=${depth[$i]} $azmthfile| gmt xyz2grd -R$R_tomo  -I0.25/0.25 -V -Gtomo_grd
#  makecpt -Z contiunous -I:inverse -M: overrule background -Di Select the back- and foreground colors, (i: from input)
             gmt makecpt -C$tomo_cpt -T$R_T -I -V >  $cpt_in
# grdimage -Q z=NaN transparent
             gmt grdimage tomo_grd -J$J -R$R_tomo -C$cpt_in -V -Bxa1f1 -Bya1f1 -BWeSn+t"$title" ${offset[$i]} -K  -O >> $pic_name

             awk '{if($3==depth1) print $1,$2,$5,$6*40 }' depth1=${depth[$i]} $azmthfile | gmt  psxy -J$J  -R$R_tomo -SV0.2c+jc  -W1.3p,black -Gred  -K -O  -V>> $pic_name

  # Add lengend
  # 0.01(1%)---*40---0.1
             # echo "101.69 23.59 1.5 0.7" | gmt psxy  -J$J -R$R_tomo -Sr0.4  -W0.8p,black  -Gwhite -K -O -V >>$pic_name
             # echo "101.45 23.59 90 0.4"  |gmt  psxy -J$J  -R$R_tomo -SV0.2c+jc  -W1.3p,black -Gred  -K -O  -V>> $pic_name
             # color=black
             # echo "101.8 23.59 0 1%" | gmt pstext -R$R_tomo -J$J -F+f20p,Times--Roman,$color+a+jCM -N  -K -O -V >> $pic_name

# 101.4/104.4/
             echo "104.05 23.65 1.5 0.7" | gmt psxy  -J$J -R$R_tomo -Sr0.4  -W0.8p,black  -Gwhite -K -O -V >>$pic_name
             echo "103.85 23.65 90 0.4"  |gmt  psxy -J$J  -R$R_tomo -SV0.2c+jc  -W1.3p,black -Gred  -K -O  -V>> $pic_name
             color=black
             echo "104.2 23.65 0 1%" | gmt pstext -R$R_tomo -J$J -F+f20p,Times--Roman,$color+a+jCM -N  -K -O -V >> $pic_name


-E sidebar triangles for back- and/or foreground colors; -I: light
             if  [ $i = 2 ]  ;then
                gmt psscale -D3.25c/-1.6c/7c/0.4ch -C$cpt_in -B$psscale_B -B"$psscale_By" -E  -K -O -V>> $pic_name
                # gmt psscale -R -J  -DjBC+w6c/0.5c+o0c/-1.8c+h+e+ml -C$cpt_in -B$psscale_B -B"$psscale_By" -K -O -V>> $pic_name
             fi
#     pan=0.5p
#     fault_color=0/0/0
#     for ((j=1; j<=$FTIndex;j=j+1));do
#        gmt psxy ${fault[j]} -J$J -R$R_tomo -W$pan,$fault_color -K -O -V >>  $pic_name
#     done
        done
        gmt psxy  -J$J -R$R_tomo -T -O -V >> $pic_name
        gmt ps2raster -P -Tf -V $pic_name
  done

  echo "finish plot."

rm tomo_grd $cpt_in gmt.history gmt.conf
rm *.ps
