if [ $# -gt 0 ]; then
    echo "  "
    echo "  ../$(basename $0)"
    echo "  "
    exit -1
fi

NJETS=0

LUMINOSITY=19365

CHANNEL="OF"

INPUTFILE="../aMCatNLO/WW0j1j_fxfx_2l2n_tau_pythia8.root"

OUTPUTFILE="../rootfiles/nominals/0jet/WWGEN/WW_GEN_0jet_amcatnlo_full.root"

root -l -b -q "prepareAMCatNLO.C(\"$INPUTFILE\",\"$OUTPUTFILE\",$NJETS,1)"
	
SYSTEMATICS="       \
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

for SYSTEMATIC in $SYSTEMATICS; do
    cp $OUTPUTFILE ../rootfiles/$SYSTEMATIC/0jet/WWGEN/WW_GEN_0jet_amcatnlo_full.root
done