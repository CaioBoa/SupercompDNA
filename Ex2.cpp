#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>

// Função para transcrever DNA em RNA
std::string transcribe_dna_to_rna(const std::string& dna_segment) {
    std::string rna_segment(dna_segment.size(), ' ');

    #pragma omp parallel for
    for (size_t i = 0; i < dna_segment.size(); i++) {
        switch (dna_segment[i]) {
            case 'A':
            case 'a':
                rna_segment[i] = 'U';
                break;
            case 'T':
            case 't':
                rna_segment[i] = 'A';
                break;
            case 'C':
            case 'c':
                rna_segment[i] = 'G';
                break;
            case 'G':
            case 'g':
                rna_segment[i] = 'C';
                break;
            default:
                rna_segment[i] = dna_segment[i];  // Caso especial, mantém o caractere
                break;
        }
    }

    return rna_segment;
}

// Função para ler o arquivo FASTA
std::string read_fasta(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line, sequence;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            if (line[0] != '>') {  // Ignora linhas de cabeçalho
                sequence += line;
            }
        }
        file.close();
    }

    return sequence;
}

// Função para salvar a sequência transcrita em um arquivo
void write_rna_to_file(const std::string& file_path, const std::string& rna_sequence) {
    std::ofstream file(file_path);
    if (file.is_open()) {
        file << rna_sequence;
        file.close();
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();  // Início da medição do tempo

    // Lista de arquivos DNA para processar
    std::vector<std::string> file_paths = {
        "fasta_files/chr1.subst.fa",
        "fasta_files/chr2.subst.fa",
        "fasta_files/chr3.subst.fa",
        "fasta_files/chr4.subst.fa",
        "fasta_files/chr5.subst.fa",
        "fasta_files/chr6.subst.fa",
        "fasta_files/chr7.subst.fa",
        "fasta_files/chr8.subst.fa",
        "fasta_files/chr9.subst.fa",
        "fasta_files/chr10.subst.fa",
        "fasta_files/chr11.subst.fa",
        "fasta_files/chr12.subst.fa",
        "fasta_files/chr13.subst.fa",
        "fasta_files/chr14.subst.fa",
        "fasta_files/chr15.subst.fa",
        "fasta_files/chr16.subst.fa",
        "fasta_files/chr17.subst.fa",
        "fasta_files/chr18.subst.fa",
        "fasta_files/chr19.subst.fa",
        "fasta_files/chr20.subst.fa",
        "fasta_files/chr21.subst.fa",
        "fasta_files/chr22.subst.fa"
        // Adicione mais arquivos conforme necessário
    };

    int files_per_process = file_paths.size() / size;
    int start_index = rank * files_per_process;
    int end_index = (rank == size - 1) ? file_paths.size() : start_index + files_per_process;

    for (int i = start_index; i < end_index; i++) {
        std::string dna_sequence = read_fasta(file_paths[i]);
        std::string rna_sequence = transcribe_dna_to_rna(dna_sequence);

        // Salvar o RNA convertido em arquivos separados
        std::string output_file = "rna_files/chr" + std::to_string(i + 1) + ".rna.fa";
        write_rna_to_file(output_file, rna_sequence);
    }

    MPI_Barrier(MPI_COMM_WORLD);  // Sincronização entre os processos
    double end_time = MPI_Wtime();  // Fim da medição do tempo

    if (rank == 0) {
        std::cout << "Tempo total de execução: " << (end_time - start_time) << " segundos." << std::endl;
    }

    MPI_Finalize();
    return 0;
}

