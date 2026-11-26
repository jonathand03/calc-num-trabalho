# Trabalho de Cálculo Numérico 2026.2

**Runge-Kutta na Simulação de Circuitos RLC:** da história do cálculo numérico à implementação
computacional. IFCE, Departamento de Telemática, prof. Francisco J. A. de Aquino.

Entrega: **1 PDF** pelo Classroom, até **04/12/2026**. Roteiro e livro da disciplina em `referencia/`.

## Estrutura

| Caminho | O que é |
|---|---|
| `trabalho_calculo_numerico.tex` | Arquivo principal. Só junta as seções. |
| `trabalho_calculo_numerico.pdf` | PDF final, versionado. Recompile antes de commitar mudança no texto. |
| `config/preambulo.tex` | Formatação (fonte, margens, legendas, referências ABNT, ambiente `codigo`). |
| `secoes/00-capa.tex`, `00-resumo.tex` | Capa, resumo, palavras-chave e sumário. |
| `secoes/01-introducao.tex` a `07-conclusao.tex` | Uma seção por arquivo, na ordem do roteiro. |
| `secoes/06-resultados/` | Uma subseção por arquivo (convergência, regimes, estabilidade, senoide, custo, `ode`, limitações). |
| `secoes/08-nota-ia.tex` | Declaração de uso de IA, exigida pelo roteiro. |
| `secoes/apendice-codigo.tex` | Apêndice A: lê o código direto dos arquivos de `codigo/`. |
| `referencias.bib` | As 5 referências. A lista ABNT sai automática. |
| `codigo/main.sce` | Entrada de dados. Roda todos os ensaios e gera resultados e figuras. |
| `codigo/rlc_*.sci` | Funções: modelo e métodos, soluções exatas, análise, ensaios e figuras. |
| `codigo/resultados/metricas_scilab.json` | Números das tabelas, gerados pelo `main.sce`. |
| `figuras/` | Gráficos em PDF, gerados pelo `main.sce`. Não editar à mão. |
| `referencia/` | Roteiro do professor e livro de Aquino (2020). Somente leitura. |

## Rodar a simulação (Scilab)

Requer o Scilab instalado (testado com a versão do Ubuntu 24.04; o código evita funções novas,
como `toJSON`, para rodar também em versões mais antigas).

**No console do Scilab**, com a pasta do trabalho como diretório atual:

```scilab
exec("codigo/main.sce", -1)
```

**No terminal**, na raiz do repositório:

```bash
scilab -nw -nb -quit -f codigo/main.sce     # gera resultados e figuras
scilab-cli -nb -quit -f codigo/main.sce     # só resultados (sem interface gráfica, sem figuras)
```

A execução regrava `codigo/resultados/metricas_scilab.json` e os 6 PDFs de `figuras/`, e mostra
um resumo no console (limites de estabilidade, ordens e comparação com o `ode`). Depois de rodar,
confira se os números das tabelas em `secoes/06-resultados/` continuam batendo com o JSON.

## Compilar o relatório

**No computador** (TeX Live com `biblatex-abnt`, `babel-portuges` e `biber`), na raiz:

```bash
latexmk
```

O `.latexmkrc` já aponta para `trabalho_calculo_numerico.tex` e roda o `biber`. Os auxiliares
(`.aux`, `.log`, `.toc`...) ficam na raiz e são ignorados pelo git. Para limpar: `latexmk -c`.

**No Overleaf:** compilador **pdfLaTeX**. Envie o repositório compactado em
**New Project → Upload Project**.

## Como citar e inserir figura

```latex
\cite{aquino2020}                 % (Aquino, 2020)
\cite[seção~10.6]{aquino2020}     % (Aquino, 2020, seção 10.6)
\textcite[p.~182]{aquino2020}     % Aquino (2020, p. 182)
```

Cite página só se ela foi conferida no livro. O que falta conferir fica marcado com
`\confirmar{...}` (aparece em vermelho no PDF). Antes da entrega, nenhuma marca pode sobrar.

Figura com legenda acima e fonte abaixo:

```latex
\begin{figure}[htbp]
\centering
\caption{Legenda da figura.}
\label{fig:exemplo}
\includegraphics[width=\textwidth]{exemplo.pdf}
\fonte{Elaborado pelos autores (2026).}
\end{figure}
```

Trecho de código no texto, no formato do livro:

```latex
\begin{codigo}
y1 = y + h * rlc_campo(t, y, p);
\end{codigo}
```

## Commits

Conventional Commits em português, com escopo. Exemplo:
`docs(resultados): discute estabilidade com Q = 50`.
