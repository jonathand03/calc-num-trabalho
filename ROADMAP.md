# Roadmap

Enviar até **03/12** (prazo final 04/12).

## Itens do roteiro do professor

**Estrutura (§4)**

- [x] [Capa](secoes/00-capa.tex#L8)
- [x] [Resumo (150 a 300 palavras) e 3 a 5 palavras-chave](secoes/00-resumo.tex#L5)
- [x] [Introdução](secoes/01-introducao.tex#L4)
- [ ] [Referencial histórico](secoes/02-referencial-historico.tex#L7): falta conferir as fontes
- [x] [Fundamentação teórica](secoes/03-fundamentacao.tex#L6)
- [x] [Algoritmos e implementação](secoes/04-algoritmos.tex#L5)
- [x] [Aplicações](secoes/05-aplicacoes.tex#L5)
- [x] [Análise dos resultados](secoes/06-resultados.tex#L6)
- [x] [Conclusão](secoes/07-conclusao.tex#L4)
- [x] [Referências ABNT](referencias.bib#L13)
- [x] [Apêndice com o código](secoes/apendice-codigo.tex#L7)

**Parte histórica (§5)**

- [x] [Tabelas, logaritmos e interpolação](secoes/02-referencial-historico.tex#L10)
- [x] [Newton, Euler, Lagrange, Gauss, Jacobi, Runge e Kutta](secoes/02-referencial-historico.tex#L16)
- [ ] [Ábaco, régua, Pascal, Leibniz e Babbage](secoes/02-referencial-historico.tex#L24): falta fonte
- [ ] [Lovelace e Babbage](secoes/02-referencial-historico.tex#L26): falta conferir
- [ ] [ENIAC e von Neumann](secoes/02-referencial-historico.tex#L28): falta conferir
- [x] [Ponto flutuante e hardware](secoes/02-referencial-historico.tex#L32)
- [x] [Linguagens e bibliotecas](secoes/02-referencial-historico.tex#L40)
- [x] [Importância atual](secoes/02-referencial-historico.tex#L46)

**Atividade computacional (§8)**

- [x] [Entrada de dados](codigo/main.sce#L27)
- [x] Método implementado: [Euler](codigo/rlc_metodos.sci#L59), [Heun](codigo/rlc_metodos.sci#L66), [RK4](codigo/rlc_metodos.sci#L77)
- [x] [Critério de parada](codigo/rlc_metodos.sci#L92)
- [x] [Cálculo do erro](codigo/rlc_analise.sci#L34)
- [x] [Tabelas](codigo/main.sce#L86) e [gráficos](codigo/main.sce#L92)
- [x] Comparação com a [solução exata](codigo/rlc_exatas.sci#L9) e com [outro método](codigo/rlc_experimentos.sci#L151)
- [x] [Precisão](secoes/03-fundamentacao.tex#L155), [convergência](codigo/rlc_experimentos.sci#L35), [estabilidade](codigo/rlc_experimentos.sci#L74) e [custo](codigo/rlc_experimentos.sci#L129)
- [ ] Testes com condição inicial diferente

**Formato (§9)**

- [x] [A4, Times 12, espaçamento 1,5, margens 2,5 cm](config/preambulo.tex#L27)
- [x] [10 a 30 páginas (tem 19)](trabalho_calculo_numerico.pdf)
- [ ] Todas as equações citadas no texto: faltam 3 ([3.2](secoes/03-fundamentacao.tex#L28), [3.4](secoes/03-fundamentacao.tex#L46), [3.7](secoes/03-fundamentacao.tex#L89))
- [x] [Código completo no apêndice](secoes/apendice-codigo.tex#L7)

**Regras (§12)**

- [ ] [Uso de IA dentro do permitido](secoes/08-nota-ia.tex#L7): decidir
- [ ] [Autores entendem o código](secoes/04-algoritmos.tex#L49)

## Falta fazer

- [ ] Citar as 3 equações soltas
- [ ] Teste com condição inicial diferente
- [ ] Conferir as 12 marcas vermelhas do PDF nos livros
- [ ] Decidir o texto da nota de IA
- [ ] Rodar o código no Scilab e saber explicar
- [ ] Revisão final e PDF
- [ ] Enviar no Classroom

## Risco e nota estimada

| Critério                         | Peso | Risco     | Nota    |
| -------------------------------- | ---- | --------- | ------- |
| Adequação ao tema                | 10%  | Baixo     | 9,5     |
| Pesquisa histórica e referências | 15%  | **Médio** | 7,0     |
| Fundamentação teórica            | 20%  | Baixo     | 9,0     |
| Implementação                    | 20%  | Baixo     | 9,0     |
| Análise dos resultados e erros   | 15%  | Baixo     | 9,0     |
| Relação com a engenharia         | 10%  | Baixo     | 8,5     |
| Organização e normas             | 10%  | Baixo     | 8,5     |
| **Total**                        |      |           | **8,7** |

Estimativa do Claude, não do professor. A conta não inclui uma possível penalização pela regra de IA.
