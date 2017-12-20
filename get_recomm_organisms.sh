#Prend un fichier assembly et un fichier avec des noms d'espece en argument et telecharge les genomes complets.
if [[ $# -eq 0 ]] ; then
    echo 'No arguments given. Needs an assembly file first and organisms file next.'
    exit 1
fi

while IFS='' read -r line2 || [[ -n "$line2" ]]; do
grep "$line2" $1 |awk -F"\t" ' $12=="Complete Genome"{print $20}' > tmp.txt #Adresses web des dossiers des especes d'archaea a genome complet de la liste donnee ds le fichier en arg 2
sed 's/\(GCA_.*\)$/\1\/\1_genomic.fna.gz/g' tmp.txt >> tmp2.txt
#cat tmp2.txt
rm tmp.txt
done < $2
while IFS='' read -r line || [[ -n "$line" ]]; do
    wget $line
    gzip -d *.gz
    rm *.gz
done < tmp2.txt
rm tmp2.txt
