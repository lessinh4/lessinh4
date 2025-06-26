#SCRIPT PARA PRE PROCESSAMENTO DE DADOS BRUTOS MICRORNAS (SE)
#FASTP PARA RETIRADA DE SEQUÊNCIAS DE BAIXA QUALIDADE, FILTRAGEM POR TAMANHO DAS SEQUÊNCIAS ASSUMINDO A DETECÇÃO AUTOMÁTICA DO FASTP

#ENTRADA E SAÍDA DE DADOS
-i / -in1
-o/-out1

#USO SIMPLES
fastp -i in.fq -o out.fq

--stdout
-e, --average_qual 33 (Phred)  #média do score de qualidade (0,05 taxa de erro)

--length_required  18
--length_limit 30 #filtro de comprimento das reads

-U --umi_loc=read1  #processamento de UMI. quando os dados possuem alto índice de duplicação

--overrepresentation_sampling-P 100-P 1   #análise de sequência super-representada

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ALINHAMENTO COM BOWTIE
#É PRECISO MODIFICAR OS DIRETÓRIOS INPUT E OUTPUT

#!/bin/bash

INPUT_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/filtering_fastp"
REF_INDEX="/home/pablooliveira/projects/Alana/genome/Homo_sapiens.GRCh38.dna.toplevel"
OUTPUT_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/alignment_bowtie"
THREADS=4

mkdir -p "$OUTPUT_DIR"

for FASTQ in "$INPUT_DIR"/cleaned_*.fastq; do
    SAMPLE=$(basename "$FASTQ" .fastq)
    echo "Processando $SAMPLE..."
    
    bowtie -x "$REF_INDEX" \
           -n 0 -l 18 -v 2 \
           --best -k 1 \
           -p "$THREADS" \
           -S \
           "$FASTQ" \
           > "${OUTPUT_DIR}/${SAMPLE}.sam" 2> "${OUTPUT_DIR}/${SAMPLE}.log"
    
    if [ -s "${OUTPUT_DIR}/${SAMPLE}.sam" ]; then
        echo "✅ ${SAMPLE}.sam gerado com sucesso!"
    else
        echo "❌ Falha no alinhamento. Log:"
        cat "${OUTPUT_DIR}/${SAMPLE}.log"
    fi
done
