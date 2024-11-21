#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <mpi.h>
#include <omp.h>
#include <iomanip>  // Para formatação do tempo

// Mapa de códons para aminoácidos representados como números
std::map<std::string, int> codon_to_amino = {
    {"UAA", 0}, {"UAG", 0}, {"UGA", 0},  // Parada
    {"AUG", 1},  // Metionina (Início)
    {"CCA", 2}, {"UCA", 2}, {"UCG", 2}, {"UCC", 2},  // Serina
    {"CAG", 3}, {"CAA", 3},  // Glutamina
    {"ACA", 5}, {"ACC", 5}, {"ACU", 5}, {"ACG", 5},  // Treonina
    {"GUG", 4}, {"GUC", 4}, {"GUA", 4}, {"GUU", 4},  // Valina
    {"UGC", 6}, {"UGU", 6},  // Treonina
    // Outros códons...
};

// Função para encontrar todas as sequências válidas de proteínas
std::vector<std::vector<std::string>> find_protein_sequences(const std::string& rna_sequence) {
    std::vector<std::vector<std::string>> proteins;

    #pragma omp parallel
    {
        std::vector<std::vector<std::string>> local_proteins;

        #pragma omp for
        for (size_t i = 0; i < rna_sequence.size() - 2; i += 3) {
            std::string codon = rna_sequence.substr(i, 3);
            if (codon == "AUG") {  // Encontrou início de uma proteína
                std::vector<std::string> protein;
                protein.push_back(codon);  // Adicionar códon de início

                bool valid_protein = false;

                for (size_t j = i + 3; j < rna_sequence.size() - 2; j += 3) {
                    std::string next_codon = rna_sequence.substr(j, 3);
                    protein.push_back(next_codon);

                    if (codon_to_amino[next_codon] == 0) {  // Encontrou códon de parada
                        valid_protein = true;
                        break;  // Finaliza a busca por esta proteína
                    }
                }

                if (valid_protein) {  // Só salva proteínas válidas
                    #pragma omp critical
                    local_proteins.push_back(protein);
                }
            }
        }

        #pragma omp critical
        proteins.insert(proteins.end(), local_proteins.begin(), local_proteins.end());
    }

    return proteins;
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
        "rna_files/chr6.rna.fa"
        // Adicione mais arquivos conforme necessário
    };

    int files_per_process = file_paths.size() / size;
    int start_index = rank * files_per_process;
    int end_index = (rank == size - 1) ? file_paths.size() : start_index + files_per_process;

    std::vector<std::vector<std::string>> local_proteins;

    for (int i = start_index; i < end_index; i++) {
        std::string rna_sequence = read_rna(file_paths[i]);
        std::vector<std::vector<std::string>> proteins = find_protein_sequences(rna_sequence);
        local_proteins.insert(local_proteins.end(), proteins.begin(), proteins.end());
    }

    double end_time = MPI_Wtime();  // Fim da medição do tempo

    if (rank == 0) {
        std::cout << "Proteínas encontradas:\n";
        for (const auto& protein : local_proteins) {
            for (const auto& codon : protein) {
                std::cout << codon << " ";
            }
            std::cout << "\n";
        }

        std::cout << "Tempo total de execução: " 
                  << std::fixed << std::setprecision(4) 
                  << (end_time - start_time) << " segundos." << std::endl;
    }

    MPI_Finalize();
    return 0;
}






