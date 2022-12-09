---
title: Curso Metagenomica Funcional 2022. Tema 4
tags: NGS, shotgun, metagenomics
description: View the slide with "Slide Mode".
---

# Curso Metagenómica Funcional 2022
Centro de Investigaciones Biológicas del Noroeste S.C.
La Paz, Baja California Sur. Diciembre 2022

**Dr. Miguel Martínez Mercado** | marmigues@gmail.com | mmartinez@cicese.mx
          

# TEMA 4. Reconstrucción de genomas *de novo* y pangenomas.

## 4.0 Ajustes y comandos útiles

Ambientes de conda que se utilizan:
```bash=
conda activate sra_tools-env
conda activate metawrap-env
conda activate DRAM_v1.4
```

Para que puedan cargar los ambientes públicos hay que revisar que tengan las siguientes líneas en el archivo `.condarc` o crearlo si no existe

```bash=
cat ~/.condarc

envs_dirs:
  - /home/aescobar/shared_envs
  - /home/kvazquez/.conda/envs
```

Revisar que los ambientes sean visibles:
```bash
conda env list
# Se despliega lista de ambientes y su ubicación.
# El ambiente activo tiene un asterisco
```

Recuerden desactivar el ambiente conda antes de desloguearse del servidor y volver a cargar el ambiente que requieran una vez que abran una nueva terminal.

Otros comandos de utilidad
```bash=
# Para enviar los procesos a correr en el background:
nohup my_command &
# En el archivo nohup.out van a encontrar el log del proceso

# Comando para revisar los procesos corriendo en background
ps -aux | grep 'user'

# Comando para matar un proceso que esta corriendo en background
kill -9 PID

# Para descargar archivos del servidor, desde su maquina
scp <usuario>@200.23.162.231:/home/user/some/path/file.txt .
```

## 4.1 Reconstrucción de genomas a partir de secuencias metagenómicas

Nos basamos en el tutorial de [metaWRAP](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md)

### 4.1.1 Datos a procesar
Para descargar secuencias crudas de metagenomas directamente del NCBI se puede hacer via *SRA-tools*. Si quisieran usar sus datos (secuencias crudas) para correr los ejercicios, se consideran los siguientes criterios de selección:
- Un metagenoma de ambiente de **baja** diversidad: ambiente extremo, fermentación, nódulos de raices, etc...
- Que la muestra no tenga mas de 50 millones de lecturas
- Secuenciación pareada, tecnología Illumina

Ejemplo de la descarga de secuencias de NCBI
```bash
**NO correr estos comandos**
conda activate /home/aescobar/shared_envs/sra_tools-env
nohup fastq-dump --split-3 SRR8197393 &
conda deactivate
```

Para efectos prácticos en el curso usaremos datos de juguete, un metagenoma artificial.
```bash=
# Crear y entrar al directorio de trabajo
mkdir -p curso_2022/data
cd curso_2022/data

# Hacer enlace
ln -s /Datos/metagenome_raw_data/*.fastq .
```


### 4.1.2 Trimming y filtrado por calidad (~20 min) 
Dentro de `curso_2022` hacemos una carpeta que se llame 'MAGs_assembly' para correr todo lo reacionado con la reconstrucción de los genomas

```bash=
conda activate metawrap-env

# Donde estoy?
pwd
# /home/ebioinfo#/curso2022

mkdir MAGs_assembly && cd MAGs_assembly && mkdir READ_QC

nohup metawrap read_qc -1 ../data/mock_1.fastq -2 ../data/mock_2.fastq -t 8 -o READ_QC/mock_reads --skip-bmtagger &
```

Al revisar el archivo nohup.out debe terminar con el mensaje:
:::success
READ QC PIPELINE COMPLETE!!!
:::

:::info
Reportar las gráficas de calidad antes y despues del procesamiento de calidad
:::

### 4.1.3 Ensamble metagenómico (~30 min).

El ensamblador que usaremos es MegaHIT.
La otra opción (metaSPAdes) requiere mucha más memoria.

```bash=
# Donde estoy?
pwd
# /home/ebioinfo#/curso2022/MAGs_assembly

# Directorio para los reads limpios
mkdir CLEAN_READS && cd CLEAN_READS

# Enlace a las secuencias limpias
ln -s ../READ_QC/mock_reads/final_pure_reads_1.fastq
ln -s ../READ_QC/mock_reads/final_pure_reads_2.fastq

# Regresamos al directorio MAGs_assembly
cd ..

# Ensamblar
nohup metawrap assembly -1 CLEAN_READS/final_pure_reads_1.fastq -2 CLEAN_READS/final_pure_reads_2.fastq -m 200 -t 8 --megahit -o ASSEMBLY &
```

:::info
Reportar estadísticas del ensamble con el script /home/aescobar/scripts/stats.pl
:::

### 4.1.5 Binning y refinamiento de bins (~2 h+)

```bash=
cd /home/<usuario>/curso_2022/MAGs_assembly/CLEAN_READS

# Correr binning (~10min)
nohup metawrap binning -o INITIAL_BINNING -t 8 -a ../ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct ./*fastq &

# Refinamiento (2h~)
nohup metawrap bin_refinement -o BIN_REFINEMENT -t 8 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 85 -x 5 &
```

:::warning
Cuando termine el proceso de refinamiento, eliminar carpeta de `INITIAL_BINNING`
:::

:::info
Reportar el número de bins de buena calidad que obtuvieron y valores de contaminación y 'completeness'
:::

## 4.2 Asignación taxonómica de genomas y árboles filogenéticos

Etiquetar taxonómicamente con GTDBtk los MAGs de alta calidad y hacer un árbol filogenético.


### 4.2.1 Organizar espacio de trabajo
```bash=
cd /home/<usuario>/curso_2022

# Crear directorio TAXONOMY
mkdir TAXONOMY && cd TAXONOMY && mkdir genomes && cd genomes

# Enlace a bins de alta calidad
ln -s ../../MAGs_assembly/CLEAN_READS/BIN_REFINEMENT/metawrap_85_5_bins/*.fa .

# Enlace a genomas de referencia que corresponden a aislados
ln -s /Datos/REF_genomes/*.fa .

cd ..
```

### 4.2.2 Flujo de trabajo GTDBtk

[GTDB](https://gtdb.ecogenomic.org/) = Genome Taxonomy Database es mantenida por el Australian Centre for Ecogenomics (ACE)
El equipo de GTDB provee un conjunto de herramientas para la asignación taxonómica de genomas de bacterias y arqueas.

El ambiente conda con las herramientas que necesitamos es `gtdbtk-2.1.1`
```bash=
# Activar ambiente
conda activate gtdbtk-2.1.1

# Identificar (~20 min)
nohup gtdbtk identify --genome_dir ./genomes --out_dir gtdb_out -x fa --cpus 4 &

# Generar alineamiento (~5 min)
nohup gtdbtk align --identify_dir gtdb_out --out_dir align_out --cpus 4 &

# Clasificar (~ 1.5 h)
nohup gtdbtk classify --genome_dir ./genomes -x fa --align_dir ./align_out --out_dir classify_out --cpus 4 &

# Descomprimir el alineamiento
gzip -d align_out/align/gtdbtk.bac120.user_msa.fasta.gz

# Inferencia filogenetica (~1 min)
nohup gtdbtk infer --msa_file ./align_out/align/gtdbtk.bac120.user_msa.fasta --out_dir infer_out --cpus 4 &

conda deactivate
```

### 4.2.3 Visualización del árbol

iTOL. 
Crear una cuenta en la siguiente dirección:
https://itol.embl.de/itol_account.cgi

Aquí encuentran la documentación necesaria para hacer las anotaciones que el programa soporta. Básicamente son tablas con los datos a cargar
https://itol.embl.de/help.cgi

Descargar el arbol a su máquina
```bash
scp <usuario>@200.23.162.231:/home/<usuario>/curso_2022/TAXONOMY/infer_out/gtdbtk.unrooted.tree .
```

Preparar archivos de anotación para el árbol:
1) Etiquetar los nodos terminales con el último nivel taxonómico de acuerdo a la anotación de GTDBtk. 
```bash=
cd TAXONOMY/classify_out

# Revisar contenido de archivo
head gtdbtk.bac120.summary.tsv

# Comando de 1-linea para obtener el archivo deseado
cat gtdbtk.bac120.summary.tsv | cut -f1,2 | sed '/user_genome/d' | sed 's/;[pcofgs]__$//;s/d__Bacteria;.*;//;s/ /_/g;s/\t[gs]__/\t/' > gtdb_labels.txt

# Revisar que resultado tenga formato adecuado
head gtdb_labels.txt
# GCF_000005845	Escherichia_coli
# GCF_000006765	Pseudomonas_aeruginosa
# GCF_000006925	Escherichia_coli
# GCF_000006945	Salmonella_enterica
# GCF_000007365	Buchnera_aphidicola_O

# Descargar gtdb_labels.txt
```

2) Indicar el contenido de GC en el anillo externo en forma de heatmap o barplot 

```bash=
# Donde estoy?
cd ~/curso_2022/TAXONOMY/genomes

# Calcular contenido de GC para 1 genoma
python /home/aescobar/scripts/gc_content.py GCF_000006765.fa

# Calcular contenido de GC de varios genomas con ciclo 'for'
ls *.fa > lista.txt
for i in $(cat lista.txt); do (python /home/aescobar/scripts/gc_content.py $i); done | sed '/^Genome/d' | cut -f1,2 > gc_content.txt

# Revisar que el resultado tenga formato adecuado
head gc_content.txt
# bin.1	61.65503662648043
# bin.2	63.22454907348612
# bin.3	39.441264039721226

# Descargar gc_content.txt
```

Subir el árbol a iTOL y agregar las anotaciones.
Indicar si el nodo terminal corresponde a un MAG o a un aislado.
Agregar leyendas.
Marcar con diferente color las ramas de cada phylum. Se puede hacer a mano.

:::success
Terminamos reconstrucción de genomas *de novo*.
:::