// BIBLIOTECAS UTILIZADAS

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <filesystem>

// FUNCOES

// imprime a rede no prompt
void imprimir_rede(std::vector<std::vector<int>> rede) {
    for (int i = 0; i < rede.size(); i++) {
        for (int j = 0; j < rede.size(); j++) {
            char caractere;
            if (rede[i][j] == -1)
                caractere = '.'; // aparencia dos spins = -1
            else
                caractere = '@'; // aparencia dos spins = 1
            std::cout << caractere << ' ';
        }
        std::cout << '\n'; // notar: podemos tratar como um soh caractere (aspas simples)
    }
    std::cout << std::endl;
};

// cria a rede quadrada
std::vector<std::vector<int>> criar_rede(int lado, double m_inicial) {
    // parametros:
    // lado = numero de celulas em cada lado da rede quadrada a ser criada
    // m_inicial = magnetizacao inicial (por spin), dentro do intervalo [-1, 1]
    
    // tamanho total da rede:
    int tamanho = lado*lado;

    // criando um vetor 1d com spins negativos:
    std::vector<int> vetor_base(tamanho, -1);

    // colocando o vetor no estado de magnetizacao inicial:
    double magnetizacao_norm = (m_inicial + 1)/2; // migrando para uma magnetizacao normalizada
                                                              // onde m = 0 corresponde a m_norm = 0.5
    int magnetizacao = std::round(magnetizacao_norm * tamanho);
    std::fill(vetor_base.begin(), vetor_base.begin() + magnetizacao, 1);

    // baguncando a ordem dos elementos no vetor:
    std::random_device gerador_semente;
    int semente = gerador_semente();
    std::mt19937 gerador_aleatorios(semente);
    std::shuffle(vetor_base.begin(), vetor_base.end(), gerador_aleatorios);

    // criando a rede (matriz quadrada):
    std::vector<std::vector<int>> resultado(lado, std::vector<int>(lado)); // cria (lado) elementos iguais a 0
    for (int i = 0; i < lado; i++) {
        for (int j = 0; j < lado; j++) {
            resultado[i][j] = vetor_base[i * lado + j]; //distribui dos elementos dos vetores em linhas
        }
    }
    return resultado;
};

// calcula m do primeiro microestado
double calcular_m_inicial(std::vector<std::vector<int>> rede) {
    double lado = rede.size();
    double soma = 0;
    for (int i = 0; i < lado; i++) {
        for (int j = 0; j < lado; j++) {
            soma += rede[i][j];
        }
    }
    double resultado = soma / (lado*lado);
    return resultado;
}

// salva uma unica matriz de rede (microestado):
void salvar_rede(std::vector<std::vector<int>> rede) {
    std::ofstream arquivo("dados-1a-tentativa/rede-inicial.csv");
    for (int linha = 0; linha < rede.size(); linha++) {
        for (int coluna = 0; coluna < rede[linha].size(); coluna++) {
            arquivo << rede[linha][coluna];
            // se nao eh o ultimo elemento da linha, add virgula
            if (coluna != rede[linha].size() - 1) {
                arquivo << ",";
            }
        }
        arquivo << "\n";
    }
    arquivo.close();
}

// salva 2 tabelas com dados de cada uma unica cadeia de Markov:
// m e estado da rede a cada passo
void salvar_dados_mc(std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > resultado_de_MC) {
    // 1o arquivo: tabela de t x m
    std::ofstream arquivo1("dados-1a-tentativa/tabela-t-m.csv");
    arquivo1 << "t,m\n";
    for (int linha = 0; linha < std::get<0>(resultado_de_MC).size(); linha++) {
        arquivo1 << std::get<0>(resultado_de_MC)[linha] << ',' << std::get<1>(resultado_de_MC)[linha] <<'\n';
    }
    arquivo1.close();
    
    // 2o arquivo: matriz com estado da rede ao final
    std::ofstream arquivo2("dados-1a-tentativa/rede.csv");
    std::vector<std::vector<int>> rede = std::get<2>(resultado_de_MC);
    for (int linha = 0; linha < rede.size(); linha++) {
        for (int coluna = 0; coluna < rede[linha].size(); coluna++) {
            arquivo2 << rede[linha][coluna];
            // se nao eh o ultimo elemento da linha, add virgula
            if (coluna != rede[linha].size() - 1) {
                arquivo2 << ",";
            }
        }
        arquivo2 << "\n";
    }
    arquivo2.close();
}

// roda a cadeia de markov n vezes, retornando uma tabela de m calculada periodicamente 
std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > rodar_MC(
    std::vector<std::vector<int>> rede,
    const double h = 0.0,
    const int n_passos = 100000,
    const int periodo_registro = 1000,
    const double T = 1.8,
    const double J = 1
) {
    // recuperando o comprimento do lado da rede:
    int lado = rede.size();
    
    // criando os vetores onde estocaremos resultados periodicos:
    std::vector<int> resultado_t = {0};
    std::vector<double> resultado_m = {calcular_m_inicial(rede)};

    // loop da cadeia de Markov:
    for (int t = 1; t < n_passos; t++) {

        // selecionando o spin candidato a virar (spin j):
        std::random_device gerador_semente; // cria o objeto
        int semente = gerador_semente(); // gera 1 semente
        std::mt19937 gerador_aleatorios(semente); // 
        std::uniform_int_distribution<> uniforme(0, lado-1);
        int linha_j = uniforme(gerador_aleatorios);
        int coluna_j = uniforme(gerador_aleatorios);

        // obtendo endereco dos vizinhos:
        // (linhas:)
        int linha_cima, linha_baixo;
        if (linha_j == 0) {
            linha_cima = lado-1;
            linha_baixo = 1;
        }
        else if (linha_j == lado-1) {
            linha_cima = lado-2; 
            linha_baixo = 0;
        }
        else {
            linha_cima = linha_j-1;
            linha_baixo = linha_j+1;
        }
        // (colunas:)
        int coluna_esquerda, coluna_direita;
        if (coluna_j == 0) {
            coluna_esquerda = lado-1;
            coluna_direita = 1;
        }
        else if (coluna_j == lado-1) {
            coluna_esquerda = lado-2; 
            coluna_direita = 0;
        }
        else {
            coluna_esquerda = coluna_j-1;
            coluna_direita = coluna_j+1;
        }

        // diferenca de hamiltonianos:
        double delta_H = 2*rede[linha_j][coluna_j]*( h + J*(
            rede[linha_cima][coluna_j] + rede[linha_baixo][coluna_j]
            + rede[linha_j][coluna_esquerda] + rede[linha_j][coluna_direita]
        ));

        // aplicando aceite:
        if (delta_H <= 0) // caso trivial, A = 1, jah viramos o spin diretamente:
            rede[linha_j][coluna_j] *= -1;
        else {
            double p = std::exp(-delta_H/T);
            std::discrete_distribution<> aceite({1-p, p});
            int A = aceite(gerador_aleatorios);
            if (A)
                rede[linha_j][coluna_j] *= -1;
        }

        // registrando estado do sistema
        if (t % periodo_registro == 0) {
            resultado_t.push_back(t);
            resultado_m.push_back(calcular_m_inicial(rede));
            // metodo otimizado parece nao estar funcionando
            // resultado_m.push_back(resultado_m.back() + 2*rede[linha_j][coluna_j]/static_cast<double>(lado*lado));
            //std::cout << "estado da rede em t = " << t << '\n';
            //imprimir_rede(rede);
        } 
    }

    // retornando resultado:
    std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > resultado(resultado_t, resultado_m, rede);
    return resultado;

};

// media e desvio padrao de uma unica cadeia de Markov:
std::vector<double> resumir_MC(std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > resultado_de_MC) {
    std::vector<double> magnetizacoes = std::get<1>(resultado_de_MC);
    std::vector<double> quarto_final(magnetizacoes.begin()+(magnetizacoes.size()*3/4), magnetizacoes.end());
    double soma = std::accumulate(quarto_final.begin(), quarto_final.end(), 0.0);
    double media = soma / quarto_final.size();
    std::vector<double> m_quadrado = quarto_final;
    for (int i = 0; i < m_quadrado.size(); i++) {
        m_quadrado[i] *= m_quadrado[i];
        m_quadrado[i] -= media*media;
    }
    double sd = std::sqrt( std::accumulate(m_quadrado.begin(), m_quadrado.end(), 0.0) / m_quadrado.size() );
    std::vector<double> resultado = {media, sd};
    return resultado;
}


// PROGRAMA

int main() {
    // parametros de criacao da rede:
    int lado = 100;
    double m_inicial = 0;

    // criando a rede:
    std::vector<std::vector<int>> rede = criar_rede(lado, m_inicial);
    // salvando a primeira rede:
    salvar_rede(rede);
    //imprimir_rede(rede);


    // simulacao:
    // vamos comecar primeiro com 1 unico loop de histerese
    // iniciamos a 1a cadeia de Markov
    std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > mc = rodar_MC(rede, 0.0);
    std::vector<double> resumo = resumir_MC(mc);
    std::vector<double> medias_m = {resumo[0]}; // para estocar m media de cada cadeia de Markov
    std::vector<double> sds_m = {resumo[1]}; // para estocar desvios padrao de m de cada cadeia de Markov
    std::vector<double> lista_h = {0.0}; // para registrar o respectivo h
    std::vector<int> lista_pernas = {1}; // para registrar a respectiva perna do loop
    // vou registrar os dados da 1a cadeia de Markov para plot
    salvar_dados_mc(mc);

    int perna = 1;
    std::cout << "p |  h  |   <m>    | sd_m\n";
    std::cout <<
            perna << " | " <<
            0.0 << " | " <<
            resumo[0] << " | " <<
            resumo[1] << '\n';
    for (double h = 0.5; h < 10; h += 0.5) {
        rede = std::get<2>(mc);
        std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > mc = rodar_MC(rede, h);
        std::vector<double> resumo = resumir_MC(mc);
        medias_m.push_back(resumo[0]);
        sds_m.push_back(resumo[1]);
        lista_h.push_back(h);
        lista_pernas.push_back(perna);
        std::cout <<
            perna << " | " <<
            h << " | " <<
            resumo[0] << " | " <<
            resumo[1] << '\n';
        // otimizando o loop para passar menos tempo na regiao de saturacao:
        if (resumo[0] > 0.95)
            break; 
    }
    perna++;
    for (double h = lista_h.back(); h > -10; h -= 0.5) {
        rede = std::get<2>(mc);
        std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > mc = rodar_MC(rede, h);
        std::vector<double> resumo = resumir_MC(mc);
        medias_m.push_back(resumo[0]);
        sds_m.push_back(resumo[1]);
        lista_h.push_back(h);
        lista_pernas.push_back(perna);
        std::cout <<
            perna << " | " <<
            h << " | " <<
            resumo[0] << " | " <<
            resumo[1] << '\n';
        if (resumo[0] < -0.95)
            break;
    }
    perna++;
    for (double h = lista_h.back(); h < 10; h += 0.5) {
        rede = std::get<2>(mc);
        std::tuple< std::vector<int>, std::vector<double>, std::vector<std::vector<int>> > mc = rodar_MC(rede, h);
        std::vector<double> resumo = resumir_MC(mc);
        medias_m.push_back(resumo[0]);
        sds_m.push_back(resumo[1]);
        lista_h.push_back(h);
        lista_pernas.push_back(perna);
        std::cout <<
            perna << " | " <<
            h << " | " <<
            resumo[0] << " | " <<
            resumo[1] << '\n';
        if (resumo[0] > 0.95)
            break;
    }     
    
} 
