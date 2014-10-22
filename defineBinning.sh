if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./$(basename $0) NJETS"
    echo "  "
    exit -1
fi

NJETS=$1


FILEPATH="root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshww/amassiro/RunI/trees/tree_skim_wwmin_09Jan2014/";

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

SAMPLES="           \
WW                  \
"


for SYSTEMATIC in $SYSTEMATICS; do

    for SAMPLE in $SAMPLES; do

	root -l -b -q "defineBinning.C(\"$FILEPATH\",$NJETS,\"$SAMPLE\",\"$SYSTEMATIC\",1)"

    done

done