# Overall
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"BurkinaFaso_chr"$2".txt"; print>"BurkinaFaso_chr"$2".txt"}' BurkinaFaso.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Cameroon_chr"$2".txt"; print>"Cameroon_chr"$2".txt"}' Cameroon.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Gambia_chr"$2".txt"; print>"Gambia_chr"$2".txt"}' Gambia.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Ghana_chr"$2".txt"; print>"Ghana_chr"$2".txt"}' Ghana.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Kenya_chr"$2".txt"; print>"Kenya_chr"$2".txt"}' Kenya.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Malawi_chr"$2".txt"; print>"Malawi_chr"$2".txt"}' Malawi.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Mali_chr"$2".txt"; print>"Mali_chr"$2".txt"}' Mali.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Nigeria_chr"$2".txt"; print>"Nigeria_chr"$2".txt"}' Nigeria.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"Tanzania_chr"$2".txt"; print>"Tanzania_chr"$2".txt"}' Tanzania.txt

# CM
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_BurkinaFaso_chr"$2".txt"; print>"CM_BurkinaFaso_chr"$2".txt"}' CM_BurkinaFaso.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Cameroon_chr"$2".txt"; print>"CM_Cameroon_chr"$2".txt"}' CM_Cameroon.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Gambia_chr"$2".txt"; print>"CM_Gambia_chr"$2".txt"}' CM_Gambia.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Ghana_chr"$2".txt"; print>"CM_Ghana_chr"$2".txt"}' CM_Ghana.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Kenya_chr"$2".txt"; print>"CM_Kenya_chr"$2".txt"}' CM_Kenya.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Malawi_chr"$2".txt"; print>"CM_Malawi_chr"$2".txt"}' CM_Malawi.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Mali_chr"$2".txt"; print>"CM_Mali_chr"$2".txt"}' CM_Mali.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Nigeria_chr"$2".txt"; print>"CM_Nigeria_chr"$2".txt"}' CM_Nigeria.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"CM_Tanzania_chr"$2".txt"; print>"CM_Tanzania_chr"$2".txt"}' CM_Tanzania.txt

# SMA
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_BurkinaFaso_chr"$2".txt"; print>"SMA_BurkinaFaso_chr"$2".txt"}' SMA_BurkinaFaso.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Cameroon_chr"$2".txt"; print>"SMA_Cameroon_chr"$2".txt"}' SMA_Cameroon.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Gambia_chr"$2".txt"; print>"SMA_Gambia_chr"$2".txt"}' SMA_Gambia.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Ghana_chr"$2".txt"; print>"SMA_Ghana_chr"$2".txt"}' SMA_Ghana.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Kenya_chr"$2".txt"; print>"SMA_Kenya_chr"$2".txt"}' SMA_Kenya.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Malawi_chr"$2".txt"; print>"SMA_Malawi_chr"$2".txt"}' SMA_Malawi.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Mali_chr"$2".txt"; print>"SMA_Mali_chr"$2".txt"}' SMA_Mali.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Nigeria_chr"$2".txt"; print>"SMA_Nigeria_chr"$2".txt"}' SMA_Nigeria.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"SMA_Tanzania_chr"$2".txt"; print>"SMA_Tanzania_chr"$2".txt"}' SMA_Tanzania.txt

# OTHER
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_BurkinaFaso_chr"$2".txt"; print>"OTHER_BurkinaFaso_chr"$2".txt"}' OTHER_BurkinaFaso.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Cameroon_chr"$2".txt"; print>"OTHER_Cameroon_chr"$2".txt"}' OTHER_Cameroon.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Gambia_chr"$2".txt"; print>"OTHER_Gambia_chr"$2".txt"}' OTHER_Gambia.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Ghana_chr"$2".txt"; print>"OTHER_Ghana_chr"$2".txt"}' OTHER_Ghana.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Kenya_chr"$2".txt"; print>"OTHER_Kenya_chr"$2".txt"}' OTHER_Kenya.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Malawi_chr"$2".txt"; print>"OTHER_Malawi_chr"$2".txt"}' OTHER_Malawi.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Mali_chr"$2".txt"; print>"OTHER_Mali_chr"$2".txt"}' OTHER_Mali.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Nigeria_chr"$2".txt"; print>"OTHER_Nigeria_chr"$2".txt"}' OTHER_Nigeria.txt
awk -F'\t' 'FNR==1{hdr=$0;next} {if (!seen[$2]++) print hdr>"OTHER_Tanzania_chr"$2".txt"; print>"OTHER_Tanzania_chr"$2".txt"}' OTHER_Tanzania.txt
