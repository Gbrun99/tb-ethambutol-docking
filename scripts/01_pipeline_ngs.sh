#!/usr/bin/env bash
set -e # Detener el script si ocurre algún error

# Definición de variables y rutas
REF="data/reference/H37Rv.fa"
R1="data/fastq/sample_R1.fastq.gz"
R2="data/fastq/sample_R2.fastq.gz"
SAMPLE="sample_tb"
THREADS=4

echo "=== 1. Control de Calidad ==="
fastqc data/fastq/*.fastq.gz -o results/

echo "=== 2. Alineamiento al Genoma de Referencia (H37Rv) ==="
# Indexar la referencia si no está indexada
[ ! -f "${REF}.bwt" ] && bwa index $REF

# Alineamiento y conversión directa a BAM ordenado
bwa mem -t $THREADS $REF $R1 $R2 | \
  samtools view -u - | \
  samtools sort -@ $THREADS -o results/bam/${SAMPLE}.sorted.bam

# Indexar archivo BAM
samtools index results/bam/${SAMPLE}.sorted.bam

echo "=== 3. Llamado de Variantes (Variant Calling) ==="
# Generar llamadas de variantes raw en VCF
bcftools mpileup -f $REF results/bam/${SAMPLE}.sorted.bam | \
  bcftools call -mv -Ob -o results/vcf/${SAMPLE}.bcf

# Filtrar por calidad de variante (QUAL > 30)
bcftools view -i 'QUAL>30' results/vcf/${SAMPLE}.bcf > results/vcf/${SAMPLE}.filtered.vcf

echo "=== Pipeline completado con éxito ==="
