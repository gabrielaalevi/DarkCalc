#!/bin/sh

homeDIR="$( pwd )"

echo "Installation will take place in $homeDIR"

cd $homeDIR


madgraph="./MadDM/MG5_aMC_v3_5_4_maddm.tar.gz"
URL=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/$madgraph
echo -n "Install MadGraph (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir MG5;
	echo "[installer] installing MadGraph5 from $madgraph"; tar -zxf $madgraph -C MG5 --strip-components 1;
	cd $homeDIR
	sed  "s|homeDIR|$homeDIR|g" mg5_configuration.txt > ./MG5/input/mg5_configuration.txt;
	echo "[installer] replacing MadDM files"
	cp MadDM/get_taacs.f MadDM/maddm.f MadDM/makefile MG5/PLUGIN/maddm/Templates/src/;
fi

cd $currentDIR
