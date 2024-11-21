#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>

// Função para contar o número de códons "AUG" (início de proteína) em uma sequência de RNA
int count_start_codons(const std::string& rna_sequence) {
    int count = 0;

    #pragma omp parallel for reduction(+:count)
    for (size_t i = 0; i < rna_sequence.size() - 2; i++) {
        if (rna_sequence[i] == 'A' && rna_sequence[i+1] == 'U' && rna_sequence[i+2] == 'G') {
            count++;
        }
    }

    return count;
}

// Função para ler o arquivo RNA
std::string read_rna(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line, sequence;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            sequence += line;
        }
        file.close();
    }

    return sequence;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();  // Início da medição do tempo

    // Lista de arquivos RNA
    std::vector<std::string> file_paths = {
        "rna_files/chr1.rna.fa",
        "rna_files/chr2.rna.fa",
        "rna_files/chr3.rna.fa",
        "rna_files/chr4.rna.fa",
        "rna_files/chr5.rna.fa",
        "rna_files/chr6.rna.fa",
        "rna_files/chr7.rna.fa",
        "rna_files/chr8.rna.fa",
        "rna_files/chr9.rna.fa",
        "rna_files/chr10.rna.fa",
        "rna_files/chr11.rna.fa",
        "rna_files/chr12.rna.fa",
        "rna_files/chr13.rna.fa",
        "rna_files/chr14.rna.fa",
        "rna_files/chr15.rna.fa",
        "rna_files/chr16.rna.fa",
        "rna_files/chr17.rna.fa",
        "rna_files/chr18.rna.fa",
        "rna_files/chr19.rna.fa",
        "rna_files/chr20.rna.fa",
        "rna_files/chr21.rna.fa",
        "rna_files/chr22.rna.fa"
        // Adicione mais arquivos conforme necessário
    };

    int files_per_process = file_paths.size() / size;
    int start_index = rank * files_per_process;
    int end_index = (rank == size - 1) ? file_paths.size() : start_index + files_per_process;

    int local_count = 0;

    for (int i = start_index; i < end_index; i++) {
        std::string rna_sequence = read_rna(file_paths[i]);
        local_count += count_start_codons(rna_sequence);
    }

    int global_count = 0;

    // Reduzir as contagens locais para o processo mestre
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);  // Sincronização antes de medir o tempo total
    double end_time = MPI_Wtime();  // Fim da medição do tempo

    if (rank == 0) {
        std::cout << "Número total de proteínas inicializadas (AUG): " << global_count << std::endl;
        std::cout << "Tempo total de execução: " << (end_time - start_time) << " segundos." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
