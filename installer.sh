#!/bin/sh

homeDIR="$( pwd )"

echo "Installation will take place in $homeDIR"

cd $homeDIR


madgraph="MG5_aMC_v3.5.6.tar.gz"
URL=https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/$madgraph
echo -n "Install MadGraph (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	mkdir MG5;
	echo "[installer] getting MadGraph5 ($URL)"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $madgraph -C MG5 --strip-components 1;
	cd $homeDIR
	cd ./MG5/bin;
	echo "[installer] installing MadDM under MadGraph5"
        echo "install maddm\nexit\n" > mad_install.txt;
	./mg5_aMC -f mad_install.txt
	cd $homeDIR
	sed  "s|homeDIR|$homeDIR|g" mg5_configuration.txt > ./MG5/input/mg5_configuration.txt;
	rm $madgraph;
	echo "[installer] replacing MadDM files"
	cp MadDM/get_taacs.f MadDM/maddm.f MadDM/makefile MG5/PLUGIN/maddm/Templates/src/;
fi

cd $currentDIR
