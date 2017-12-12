#Prend un fichier assembly en argument et telecharge les genomes complets.
echo "L'argument que vous avez entre est : "$1
awk -F"\t" ' $12=="Complete Genome"{print $20}' $1 >tmp.txt #Adresses web des dossiers des especes d'archaea a genome complet
#sed 's/[/]GCA_.*$/_genomic.fna.gz/g' tmp.txt 
sed 's/\(GCA_.*\)$/\1\/\1_genomic.fna.gz/g' tmp.txt > tmp2.txt
rm tmp.txt
while IFS='' read -r line || [[ -n "$line" ]]; do
    wget $line
    gzip -d *.gz
    rm *.gz
done < tmp2.txt
rm tmp2.txt
