#SCRIPT PARA PRE PROCESSAMENTO DE DADOS BRUTOS MICRORNAS (SE)
__________________________________________________________________________________________________________________________________________________________________________________________
#FASTP PARA RETIRADA DE SEQUÊNCIAS DE BAIXA QUALIDADE, FILTRAGEM POR TAMANHO DAS SEQUÊNCIAS ASSUMINDO A DETECÇÃO AUTOMÁTICA DO FASTP
#!/bin/bash

# Define os diretórios de entrada, saída e relatórios (PODEM SER MODIFICADOS)
INPUT_DIR="/home/pablooliveira/projects/Alana/GSE166160/"
OUTPUT_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/filtrados"
REPORTS_DIR="/home/pablooliveira/projects/Alana/GSE166160/preprocessing_data/filtrados/relatorios"

# Define o número de threads (núcleos de CPU) para o fastp
THREADS=4

# --- Início da Execução do Script ---

echo "Iniciando o script de pré-processamento de FASTQ com fastp (saída descompactada)..."

# Cria os diretórios de saída e relatórios se eles não existirem
echo "Criando diretórios de saída se não existirem..."
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORTS_DIR"
echo "Diretórios criados: $OUTPUT_DIR e $REPORTS_DIR"
echo ""

# Habilita o nullglob para evitar erros se não houver arquivos correspondentes
shopt -s nullglob
fastq_files=("$INPUT_DIR"/*.fastq* "$INPUT_DIR"/*.fq*)
shopt -u nullglob

# Verifica se algum arquivo FASTQ foi encontrado
if [ ${#fastq_files[@]} -eq 0 ]; then
    echo "❌ Erro: Nenhum arquivo FASTQ encontrado em '$INPUT_DIR'."
    exit 1
fi

echo "Iniciando o pré-processamento dos arquivos single-end com fastp..."
echo "Configurações de fastp:"
echo "  - Saída descompactada (.fastq)"
echo "  - Qualidade mínima (Phred): 33"
echo "  - Comprimento da read: 18-30bp (ideal para miRNAs)"
echo "  - Threads: $THREADS"
echo ""

# Loop através de cada arquivo FASTQ
for fastq_file in "${fastq_files[@]}"
do
    echo "--- Processando arquivo: $(basename "$fastq_file") ---"

    # Extrai o nome base do arquivo (sem extensões)
    base_name=$(basename "$fastq_file")
    base_name="${base_name%.gz}"
    base_name="${base_name%.fastq}"
    base_name="${base_name%.fq}"

    # Comando fastp:
   fastp -i "$fastq_file" -o "$OUTPUT_DIR/${base_name}_cleaned.fastq" -h "$REPORTS_DIR/${base_name}.html" -j "$REPORTS_DIR/${base_name}.json" -e 33 --length_required 18 --length_limit 30 -p -w "$THREADS"

    # Verifica o status de saída do fastp
    if [ $? -eq 0 ]; then
        echo "✅ Pré-processamento de '${base_name}' concluído! Arquivo descompactado gerado: ${base_name}_cleaned.fastq"
    else
        echo "❌ Erro no processamento de '${base_name}'."
    fi
    echo ""
done

echo "🎉 Todos os arquivos foram processados."
echo "Arquivos descompactados em: $OUTPUT_DIR"
echo "Relatórios em: $REPORTS_DIR"

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
