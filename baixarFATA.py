import os
import subprocess
from urllib import request

def download_and_extract_fasta():
    base_url = "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/snp147Mask/"
    file_template = "chr{}.subst.fa.gz"
    chromosomes = list(range(1, 23))  # Cromossomos de 1 a 22

    # Cria um diretório para os arquivos FASTA, se não existir
    os.makedirs("fasta_files", exist_ok=True)

    for chrom in chromosomes:
        file_name = file_template.format(chrom)
        file_url = base_url + file_name
        local_file_path = os.path.join("fasta_files", file_name)

        # Baixa o arquivo .gz
        print(f"Baixando {file_name}...")
        request.urlretrieve(file_url, local_file_path)

        # Descompacta o arquivo
        print(f"Descompactando {file_name}...")
        subprocess.run(["gunzip", "-k", local_file_path])  # Mantém o arquivo original .gz

        print(f"{file_name} descompactado com sucesso!")

if __name__ == "__main__":
    download_and_extract_fasta()
