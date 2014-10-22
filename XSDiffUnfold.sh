if [ $# -lt 1 ]; then
    echo "  "
    echo "  ../$(basename $0) NJETS"
    echo "  "
    exit -1
fi

NJETS=$1

LUMINOSITY=19365

CHANNEL="OF"

FILEPATH="rootfiles/"

SYSTEMATICS="       \
nominals            \
JER_down            \
JER_up              \
electronResolution  \
electronScale_down  \
electronScale_up    \
jetEnergyScale_down \
jetEnergyScale_up   \
metResolution       \
metScale_down       \
metScale_up         \
muonResolution      \
muonScale_down      \
muonScale_up        \
"

DIFFERENTIALS="\
0              \
1              \
2              \
3              \
"


for SYSTEMATIC in $SYSTEMATICS; do

    for DIFFERENTIAL in $DIFFERENTIALS; do

	root -l -b -q "XSDiffUnfold.C($LUMINOSITY,$NJETS,\"$CHANNEL\",\"$FILEPATH/$SYSTEMATIC\",false,0,false,$DIFFERENTIAL)"
	
    done

done