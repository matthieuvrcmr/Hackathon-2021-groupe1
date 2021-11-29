# Hackathon-2021-groupe1

#Contacts : 

# donatien.lege@agroparistech.fr
# mathilde.guyot@agroparistech.fr
# matthieu.vercaemer@agroparistech.fr
# anthony.boutard@universite-paris-saclay.fr

# Si vous souhaitez lancer notre pipeline suivez les instructions suivantes :

# Placez vous dans une machine virtuelle 

# Se placer dans mydatalocal :
cd /mnt/mydatalocal

# Activer l'environnement Bioconda :
conda activate 

# Se connecter à notre repository Git : 
git clone git@github.com:DonatienLege/Hackathon-2021-groupe1

# Se mettre dans le dossier contenant le Workflow nommé "Snakefile" : 
cd Hackathon-2021-groupe1/Workflow

# Activer snakemake : 
conda activate snakemake

# Installer singularity : 
conda install singularity

# Exécuter le workflow avec la commande suivante : 
snakemake --use-singularity --cores all
