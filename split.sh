#Prend un fichier assembly en argument et telecharge les genomes complets.
if [[ $# -eq 0 ]] ; then
    echo 'No arguments given. Needs an assembly file.'
    exit 1
fi

awk -F"\t" ' $12=="Complete Genome"{print $20}' $1 >tmp.txt #Adresses web des dossiers des especes d'archaea a genome complet
sed 's/\(GCA_.*\)$/\1\/\1_genomic.fna.gz/g' tmp.txt > tmp2.txt
rm tmp.txt
while IFS='' read -r line || [[ -n "$line" ]]; do
    wget $line
    gzip -d *.gz
    rm *.gz
done < tmp2.txt
rm tmp2.txt
