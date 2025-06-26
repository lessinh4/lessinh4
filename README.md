#SCRIPT PARA PRE PROCESSAMENTO DE DADOS BRUTOS MICRORNAS (SE)
#FASTP PARA RETIRADA DE SEQUÊNCIAS DE BAIXA QUALIDADE, FILTRAGEM POR TAMANHO DAS SEQUÊNCIAS ASSUMINDO A DETECÇÃO AUTOMÁTICA DO FASTP

#INSTALAR FASTP BIOCONDA
conda install -c bioconda fastp

#ENTRADA E SAÍDA DE DADOS
-i / -in1
-o/-out1

#USO SIMPLES
fastp -i in.fq -o out.fq

--stdout

#média do score de qualidade (0,05 taxa de erro)
-e, --average_qual 33 (Phred)

#filtro de comprimento das reads
--length_required  18
--length_limit 30

#processamento de UMI
#quando os dados possuem alto índice de duplicação
-U --umi_loc=read1

#análise de sequência super-representada
--overrepresentation_sampling-P 100-P 1
