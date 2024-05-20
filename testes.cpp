// esse arquivo serve para testar uma serie de funcionalidades do c++ antes de comecar o projeto,
// visto que ainda estou aprendendo a linguagem

// sobre arrays vs vectors:
// ok, vou me render a vectors...
// 

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

void imprimir_vetor(std::vector<int> vetor) {
    for (int i = 0; i < vetor.size(); i++) {
        std::cout << vetor[i];
    }
    std::cout << std::endl;
};

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

std::vector<std::vector<int>> criar_rede(int lado, double magnetizacao_por_spin) {
    // tamanho total da rede:
    int tamanho = lado*lado;

    // criando um vetor 1d por enquanto de forma estruturada (nao-aleatoria),
    // que sera a base para criacao da rede:
    std::vector<int> vetor_base(tamanho, -1);

    // colocando esse vetor no estado de magnetizacao dado pelo usuario:
    double magnetizacao_norm = (magnetizacao_por_spin + 1)/2; // migrando para uma magnetizacao normalizada
                                                              // onde m = 0 corresponde a m_norm = 0.5
    int magnetizacao = std::round(magnetizacao_norm * tamanho);
    std::fill(vetor_base.begin(), vetor_base.begin() + magnetizacao, 1);
    // apenas visualizando:
    // imprimir_vetor(vetor_base);

/*
    // (parenteses: apenas testando algumas coisas para aprender:
    std::vector<int>::iterator t = vetor_base.begin(); // notar, se eu quero ver o que sai nessa begin
                                                       // preciso saber de antemao seu tipo, que eh
                                                       // esse iterador ai
    int* p = &(*t); // eh ilicito imprimir o iterador entao vamos pegar seu endereco de memoria
                    // notar que podemos desrefenciar iteradores tal como ponteiros
    // imprimindo o endereco de memoria do iterador:
    std::cout << "[teste - endereco do iterador: " << p << "]" << std::endl; 
    // imprimindo o valor parao qual o iterador aponta:
    std::cout << "[teste - valor apontado pelo iterador: " << *t << "]" << std::endl;
    // enfim, vamos voltar ao codigo )
*/

    // agora queremos baguncar a ordem dos elementos no vetor
  
    // gerador de seed baseado em configuracoes de entropia do sistema (n entendi 100%):
    std::random_device gerador_semente; // cria o objeto
    int semente = gerador_semente(); // gera 1 semente
    // apenas para entender:
    // std::cout << semente << std::endl;
    
    // criacao de um objeto gerador de numeros aleatorios, atribuindo o nome gerador_aleatorios
    // fazendo gerador_aleatorios(semente), estamos chamando o construtor desse objeto,
    // ja quando fazemos gerador_aleatorios(), estamos gerando 1 numero aleatorio
    std::mt19937 gerador_aleatorios(semente);
    // apenas para entender:
    // std::cout << gerador_aleatorios() << std::endl;
    // std::cout << gerador_aleatorios() << std::endl;

    // vamos usar shuffle da biblioteca <algorithm>
    // ele precisa de um argumento de comeco, um argumento de fim (ambos iteradores)
    // e um argumento de objeto de gerador aleatorio (este objeto ja tem todas as
    // especificacoes necessarias para shuffle gerar numeros suficientes preenchendo
    // os indices de cada elemento do vetor, pelo que entendi)
    std::shuffle(vetor_base.begin(), vetor_base.end(), gerador_aleatorios);
    // apenas visualizando:
    // imprimir_vetor(vetor_base);

    // vamos agora criar a rede (uma matriz quadrada basicamente),
    // que vai ser simplesmente um vetor de vetores
    // eh melhor que usar array ou std::array porque eh dinamico
    std::vector<std::vector<int>> resultado(lado, std::vector<int>(lado)); // cria (lado) elementos, todos com
                                                                      // valor = vetor de lado zeros
    for (int i = 0; i < lado; i++) {
        for (int j = 0; j < lado; j++) {
            resultado[i][j] = vetor_base[i * lado + j]; //distribui dos elementos dos vetores em linhas
        }
    }
    return resultado;
};

double calcular_magnetizacao_por_spin(std::vector<std::vector<int>> rede) {
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


// funcao que roda a cadeia de markov n vezes,
// retornando uma tabela com a magnetizacao 
std::pair< std::vector<int>, std::vector<double> > rodar_MC(
    std::vector<std::vector<int>> rede,
    const double m_por_spin_inicial,
    const double h = 0.0,
    const int n_passos = 1000,
    const int periodo_registro = 100,
    const double T = 2.0,
    const double J = 1.0
) {
    // apenas recuperando o comprimento do lado da rede:
    int lado = rede.size();
    
    // criando os vetores onde estocaremos resultados periodicos:
    std::vector<int> resultado_t = {0};
    std::vector<double> resultado_m = {m_por_spin_inicial};

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
            resultado_m.push_back(calcular_magnetizacao_por_spin(rede));
            //std::cout << "estado da rede em t = " << t << '\n';
            //imprimir_rede(rede);
        } 
    }

    // retornando resultado:
    std::pair< std::vector<int>, std::vector<double> > resultado = {resultado_t, resultado_m};
    return resultado;

};


int main() {
    // vamos pegar do usuario o tamanho do modelo (lado da rede quadrada):
    std::cout << "Digite o comprimento do lado da rede: ";
    int lado;
    std::cin >> lado;

    // idem a magnetizacao por spin:
    // percentual de spins = 1
    std::cout << "Digite a magnetizacao por spin do modelo: ";
    double magnetizacao_por_spin;
    std::cin >> magnetizacao_por_spin;

    // criando a rede:
    std::vector<std::vector<int>> rede = criar_rede(lado, magnetizacao_por_spin);
    //imprimir_rede(rede);

    // teste: 
    // rodando os passos de Monte Carlo para um único campo externo h
    const double h = 0;
    const double J = 0.75;
    const double T = 2;
    const int n_passos = 1000000;
    const int periodo_registro = 10000;
    std::pair< std::vector<int>, std::vector<double>> resultados = rodar_MC(rede, magnetizacao_por_spin, h, n_passos, periodo_registro, T, J);

    // tabela com os resultados:
    std::cout << "t    |     m\n";
    for (int i = 0; i < resultados.first.size(); i++) {
         std::cout << resultados.first[i] << " | "<< resultados.second[i] << '\n';
    }

    // analise dos resultados:
    std::vector<double> metade(resultados.second.begin()+(resultados.second.size()/2), resultados.second.end());
    double soma = std::accumulate(metade.begin(), metade.end(), 0.0);
    double media = soma / metade.size();
    std::vector<double> m_quadrado = metade;
    for (int i = 0; i < m_quadrado.size(); i++) {
        m_quadrado[i] *= m_quadrado[i];
        m_quadrado[i] -= media*media;
    }
    double sd = std::sqrt( std::accumulate(m_quadrado.begin(), m_quadrado.end(), 0.0) / m_quadrado.size() );
    std::cout << "<m> = " << media << " | sd = " << sd;
}