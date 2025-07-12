#SCRIPT PARA PRE PROCESSAMENTO DE DADOS BRUTOS MICRORNAS (SE)
#FASTP PARA RETIRADA DE SEQUÊNCIAS DE BAIXA QUALIDADE, FILTRAGEM POR TAMANHO DAS SEQUÊNCIAS ASSUMINDO A DETECÇÃO AUTOMÁTICA DO FASTP

#ENTRADA E SAÍDA DE DADOS
-i / -in1
-o/-out1

#USO SIMPLES
fastp -i in.fq -o out.fq


-e, --average_qual 33 (Phred)  #média do score de qualidade (0,05 taxa de erro)
--stdout

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
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#RELATÓRIOS DE QUALIDADE SAMTOOLS E QUALIMAP
#!/bin/bash
# Configurações básicas modificáveis
SAM_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/alignment_bowtie"
OUTPUT_DIR="$SAM_DIR/quality_reports"
mkdir -p "$OUTPUT_DIR"

# Loop para processar cada arquivo .sam
for SAM_FILE in "$SAM_DIR"/*.sam; do
    SAMPLE=$(basename "$SAM_FILE" .sam)
    echo "Processando $SAMPLE..."
    
    # 1. Converter SAM para BAM (formato binário)
    samtools view -Sb "$SAM_FILE" > "$OUTPUT_DIR/${SAMPLE}.bam"
    
    # 2. Ordenar BAM (necessário para estatísticas)
    samtools sort "$OUTPUT_DIR/${SAMPLE}.bam" -o "$OUTPUT_DIR/${SAMPLE}.sorted.bam"
    
    # 3. Gerar estatísticas básicas (flagstat)
    samtools flagstat "$OUTPUT_DIR/${SAMPLE}.sorted.bam" > "$OUTPUT_DIR/${SAMPLE}_flagstat.txt"
    
    echo "✅ $SAMPLE: Conversão e flagstat concluídos!"
done

echo "Análise finalizada! Verifique os arquivos em $OUTPUT_DIR"

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#GERANDO MATRIZ DE CONTAGEM FEATURE COUNTS miRNA SINGLE END
#!/bin/bash
# Ativa o ambiente Conda se necessário
# source activate bioconda_env  # Descomente se estiver fora do ambiente

# Caminhos (MODIFICAR QUANDO FOR USAR)
BAM_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/quality_reports/bam_files"
ANNOTATION="/home/pablooliveira/projects/Alana/hsa.gff3"
OUT_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/counts_matrix"
POST_MATRIX_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/quality_reports/post_matrix"

# Garante que diretórios existem
mkdir -p "$OUT_DIR"
mkdir -p "$POST_MATRIX_DIR"

# Coleta todos os arquivos *.sorted.bam
BAM_FILES=$(ls $BAM_DIR/*sorted.bam)

# Roda featureCounts
featureCounts \
  -T 8 \
  -M --primary \
  -t miRNA \
  -g Name \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/matrix_counts.txt" \
  $BAM_FILES

# Checa se houve sucesso
if [ $? -eq 0 ]; then
  echo "Contagem gerada com sucesso em $OUT_DIR/matrix_counts.txt"
else
  echo "Erro ao executar featureCounts" >&2
fi
