#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cctype>
#include <mpi.h>
#include <omp.h>

std::unordered_map<char, int> count_bases(const std::string& sequence) {
    std::unordered_map<char, int> global_counts = {{'a', 0}, {'t', 0}, {'c', 0}, {'g', 0}};
    
    // Vetores temporários para contagem paralela
    int local_counts[4] = {0, 0, 0, 0};  // Índices: 0->'a', 1->'t', 2->'c', 3->'g'

    #pragma omp parallel
    {
        int thread_counts[4] = {0, 0, 0, 0};

        #pragma omp for
        for (size_t i = 0; i < sequence.size(); i++) {
            char base = tolower(sequence[i]);
            switch (base) {
                case 'a': thread_counts[0]++; break;
                case 't': thread_counts[1]++; break;
                case 'c': thread_counts[2]++; break;
                case 'g': thread_counts[3]++; break;
                default: break;  // Ignorar caracteres não relevantes
            }
        }

        #pragma omp critical
        {
            for (int j = 0; j < 4; j++) {
                local_counts[j] += thread_counts[j];
            }
        }
    }

    // Atualizar o mapa global com os valores finais
    global_counts['a'] = local_counts[0];
    global_counts['t'] = local_counts[1];
    global_counts['c'] = local_counts[2];
    global_counts['g'] = local_counts[3];

    return global_counts;
}

std::string read_fasta(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line, sequence;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            if (line[0] != '>') {
                sequence += line;
            }
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

    double global_start_time = MPI_Wtime();

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

    std::unordered_map<char, int> global_counts = {{'a', 0}, {'t', 0}, {'c', 0}, {'g', 0}};
    std::unordered_map<char, int> local_counts;

    double local_start_time = omp_get_wtime();

    for (int i = start_index; i < end_index; i++) {
        std::string sequence = read_fasta(file_paths[i]);
        local_counts = count_bases(sequence);

        for (auto& pair : local_counts) {
            global_counts[pair.first] += pair.second;
        }
    }

    double local_end_time = omp_get_wtime();
    std::cout << "Tempo de execução local (processo " << rank << "): " 
              << (local_end_time - local_start_time) << " segundos." << std::endl;

    std::unordered_map<char, int> final_counts = {{'a', 0}, {'t', 0}, {'c', 0}, {'g', 0}};
    for (auto& pair : global_counts) {
        int total;
        MPI_Reduce(&pair.second, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            final_counts[pair.first] = total;
        }
    }

    if (rank == 0) {
        double global_end_time = MPI_Wtime();
        std::cout << "Contagem Final:\n";
        std::cout << "A: " << final_counts['a'] << "\n";
        std::cout << "T: " << final_counts['t'] << "\n";
        std::cout << "C: " << final_counts['c'] << "\n";
        std::cout << "G: " << final_counts['g'] << "\n";
        std::cout << "Tempo total de execução: " 
                  << (global_end_time - global_start_time) << " segundos." << std::endl;
    }

    MPI_Finalize();
    return 0;
}




