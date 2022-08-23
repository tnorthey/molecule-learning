
cd generated

BAGEL="/home/thomas/bagel/obj/src/BAGEL"
for i in {0000..9999}
do
	$BAGEL "$i.json" > "$i.out"
done

cd -
