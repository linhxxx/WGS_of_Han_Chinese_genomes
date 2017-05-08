#argv1 hg18 spratio
#argv2 scaf spratio

rm -f hg18.instruction scaf.instruction


echo "2" >>hg18.instruction
echo "hg18.distribution" >>hg18.instruction
echo "3" >> hg18.instruction
echo "insertion_hg18.list" >> hg18.instruction
echo "insertion_hg18.picked" >>hg18.instruction
echo "3" >>hg18.instruction
echo "deletion_hg18.list" >>hg18.instruction
echo "deletion_hg18.picked" >>hg18.instruction
echo "3" >>hg18.instruction
echo "inversion_hg18.list" >>hg18.instruction
echo "inversion_hg18.picked" >>hg18.instruction
echo "0" >>hg18.instruction

echo "2" >>scaf.instruction
echo "scaf.distribution" >>scaf.instruction
echo "3" >> scaf.instruction
echo "insertion_scaf.list" >> scaf.instruction
echo "insertion_scaf.picked" >>scaf.instruction
echo "3" >>scaf.instruction
echo "deletion_scaf.list" >>scaf.instruction
echo "deletion_scaf.picked" >>scaf.instruction
echo "3" >>scaf.instruction
echo "inversion_scaf.list" >>scaf.instruction
echo "inversion_scaf.picked" >>scaf.instruction
echo "0" >>scaf.instruction

nohup scrPATH/readspratio $1 < hg18.instruction &
nohup scrPATH/readspratio $2 < scaf.instruction &


echo "Done!"

