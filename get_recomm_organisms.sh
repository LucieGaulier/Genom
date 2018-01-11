#Prend un fichier assembly et un fichier avec des noms d'espece en argument et telecharge les genomes complets.
if [[ $# -lt 2 ]] ; then
    echo 'Missing arguments. Needs an assembly file first and organisms file next.'
    exit 1
fi

rm -f tmp2.txt

while IFS='' read -r line2 || [[ -n "$line2" ]]; do
echo "species $line2"
grep "$line2" $1 |awk -F"\t" '$12=="Complete Genome"{print $20}' > tmp.txt #Adresses web des dossiers des especes d'archaea a genome complet de la liste donnee ds le fichier en arg 2 #$20 l'adresse, $8 le nom.
sed 's/\(GCA_.*\)$/\1\/\1_genomic.fna.gz/g' tmp.txt >> tmp2.txt
cat tmp.txt
rm tmp.txt
done < $2
while IFS='' read -r line || [[ -n "$line" ]]; do
    wget -nv --show-progress $line
    gzip -d *.gz
    rm -f *.gz
done < tmp2.txt
rm tmp2.txt
