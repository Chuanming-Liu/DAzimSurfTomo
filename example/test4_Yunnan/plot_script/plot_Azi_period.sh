#!/bin/bash
# gmt 5 script of ploting 3D Vs and azimuthal anisotropy slices
# Author: Chuanming Liu

cpt_path=.
input=Period_Azm_reSmp.inv
out='Az_period.ps'

# necessary files
tomo_cpt=$cpt_path/"BlueDarkRed18.cpt"
cpt_in=./tomo_mk.cpt
fault=Fault_Yunnan.txt
vsclip=YN_clicp.txt
azmclip=Az_clicp.txt

# Plot setting
# depth
period[1]=10
period[2]=20
period[3]=30
period[4]=40
Tnum=4

# GMT 5 defaults
gmt set PS_MEDIA A3
gmt set MAP_TICK_LENGTH_PRIMARY 4p/1p 
gmt set MAP_FRAME_TYPE Plain
gmt set MAP_FRAME_WIDTH 0.11c
gmt set MAP_FRAME_PEN 1.5p
gmt set FONT Times-Roman
gmt set FONT_ANNOT_PRIMARY 21p,Times-Roman 
gmt set FONT_TITLE 22p,Times-Roman
gmt set FONT_LABEL 19p,Times-Roman  
gmt set MAP_TITLE_OFFSET 0.0c
gmt set MAP_GRID_PEN_PRIMARY thinner,black
gmt set COLOR_BACKGROUND 0/0/255
gmt set COLOR_FOREGROUND 255/0/0
gmt set COLOR_NAN 255/255/255
gmt set COLOR_MODEL RGB


offset[1]="-Y16c -X1c"
offset[2]=-X11c
offset[3]="-X-11c -Y-11c"
offset[4]="-X11c"

label=("(a)"  "(a)" "(b)" "(c)" "(d)" "(e)" "(f)" )


J=m1c
R=98/107/21/29
Range=-8/8/0.01

gmt makecpt -C$tomo_cpt -I -T$Range > tomo.cpt
gmt psxy -J$J -R$R -T  -K > $out
for ((i=1;i<=$Tnum;i=i+1));do
     title="period "${period[$i]}
     echo $title
     awk '{if($3==period1) print $1,$2,$4}' period1=${period[$i]}  $input| gmt surface -R$R -I0.01  -Gtomo.grd
     gmt psbasemap -J$J -R$R -Bxa2f2 -Bya2f2 -BWeSn ${offset[$i]} -K -O>>$out
     gmt psclip $vsclip -J -R -K -O>> $out
     gmt grdimage tomo.grd -J -R -C"tomo.cpt"  -K -O>> $out
     gmt psclip -C -K -O>> $out
     color=105/105/105
     gmt psxy $fault -J -R -W1.6p,$color -K -O>>  $out

     # plot anisotropy
     gmt psclip $azmclip -J -R -K -O>> $out
     awk '{if($3==pd1) print $1,$2,$6,$7*40 }' pd1=${period[$i]} $input | gmt  psxy -J$J  -R$R -SV0.2c+jc -W2p,black -Gred  -K -O>> $out
     gmt psclip -C -K -O  >> $out

     gmt pslegend -J -R  -DjTL+w1c/0.8c+o0c/0c+l0.9 -C0p/0p -F+gwhite+p1.1p,black -K -O<< EOF >> $out
L 19 C ${label[$i]}
EOF

    gmt pslegend -J -R -C0p/0p -DjBR+w1.8c/0.8c+o0c/0c+l0.9  -F+gwhite+p1.1p,black -K -O<< EOF >> $out
S 0.4c v0.2c+jc 0.4c black 2p,black 0.8c 1%
EOF

     word="T "${period[$i]}" s"
     gmt pslegend -J -R  -DjTR+w3.2c/0.8c+o0c/0c+l0.9 -C0p/0p -F+gwhite+p1.1p,black -K -O<< EOF >> $out
L 21 C $word
EOF

    if  [ $i = 5 ] ;then
       gmt psscale -D4c/-1.6c/7c/0.4ch -C"tomo.cpt" -Bxa0.4f0.1 -By+l"Vs (km/s)"  -K -O>> $out
    fi
done

gmt psxy  -J$J -R$R -T -O  >> $out
gmt psconvert -P -Tf -A0.5c $out

rm tomo.grd tomo.cpt gmt.history
rm gmt.conf
rm *.ps
