# simulator-for-benchmark-datasets
Bachelorthesis. Prototyp für Simulation für synthetische genomischen Benchmarkdatensätzen

# Testen
Für das Testen bitte die Rohdaten in einem Unterordner /Rohdaten bereitstellen. Benötigt werden eine VCF-Datei filteredonlychm1.vcf bzw. komprimiert und indixiert filteredonlychm1.vcf.gz.
Dafür die entsprechenden Dateien unter https://github.com/lh3/CHM-eval/releases herunterladen und mit VCFtools nur noch die Varianten für Chromosom21 herausgefiltert

vcftools --vcf input_file.vcf --chr 21 --out onlychm21

Anschließend mit der Bedtools intersect mit der BED-Datei abgleichen: 

bedtools intersect -a filteredonlychm1.vcf -b full.38.bed > filtered21.vcf

Mit dieser VCF-Datei startet man nun eine Simulation mit einer gewünschen Variantenfrequenz von bspw. 0.01, wobei man das Schreiben der BAM-Datei aus Gründen der Schnelligkeit auskommentieren kann.
Die dabei generierte VCF-Datei nennt man dann filteredonlychm1 und nachdem sie komprimiert und mittels vcftools indixiert wurde, diente sie so auch als Grundlage für die Laufzeitmessungen. Alternativ kann man auch bcftools verwenden.

bgzip -c filteredonlychm1.vcf > filteredonlychm1.vcf.gz

vcftools index filteredonlychm1.vcf.gz 


Die BAM-Datei kann mittels Samtools view heruntergeladen und auf Chr21 beschränkt werden: 

samtools view -b <PathtoBam> chr21 > inchr21.bam
  
Auf https://www.ebi.ac.uk/ena/browser/view/PRJEB13208 gibt es die verschiedenen BAM-Dateien zum Download.
Die Indixierung erfolgte mittels Samtools: 
  
samtools index inchr21.bam

# License
Published under the MIT License
