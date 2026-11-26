# Trabalho de Cálculo Numérico 2026.2

História do Cálculo Numérico e sua relação com a Ciência da Computação, com aprofundamento nos
métodos de Runge-Kutta aplicados ao circuito RLC série. IFCE, Departamento de Telemática,
prof. Francisco J. A. de Aquino.

Entrega: **1 PDF** pelo Classroom, até **04/12/2026**. Enunciado em `referencia/`.

## Estrutura

| Caminho | O que é |
|---|---|
| `trabalho_calculo_numerico.tex` | Arquivo principal. Só junta as seções. |
| `trabalho_calculo_numerico.pdf` | PDF final, versionado. Recompile antes de commitar mudança no texto. |
| `config/preambulo.tex` | Formatação (fonte, margens, legendas, referências ABNT). |
| `secoes/00-capa.tex`, `00-resumo.tex` | Capa, resumo, palavras-chave e sumário. |
| `secoes/01-introducao.tex` a `07-conclusao.tex` | Uma seção por arquivo, na ordem do enunciado. |
| `secoes/06-resultados/` | Uma subseção por arquivo (convergência, regimes, estabilidade etc.). |
| `secoes/08-nota-ia.tex` | Declaração de uso de IA, exigida pelo enunciado. |
| `secoes/apendice-codigo.tex` | Apêndice A: lê o código direto de `codigo/rlc_runge_kutta.py`. |
| `referencias.bib` | Referências. A lista ABNT sai automática. |
| `codigo/rlc_runge_kutta.py` | Implementação dos métodos, experimentos e figuras. |
| `codigo/resultados/metricas.json` | Números das tabelas, gerados pelo script. |
| `figuras/` | Gráficos em PDF, gerados pelo script. Não editar à mão. |
| `referencia/` | Enunciado do trabalho. |

## Compilar o relatório

**No computador** (TeX Live com `biblatex-abnt`, `babel-portuges` e `biber`), na raiz:

```bash
latexmk
```

O `.latexmkrc` já aponta para `trabalho_calculo_numerico.tex` e roda o `biber`. Os auxiliares
(`.aux`, `.log`, `.toc`...) ficam na raiz e são ignorados pelo git. Para limpar: `latexmk -c`.

**No Overleaf:** compilador **pdfLaTeX**. Envie o repositório compactado em
**New Project → Upload Project**.

## Gerar figuras e métricas

Requer Python 3 com `numpy` e `matplotlib`.

```bash
python3 codigo/rlc_runge_kutta.py
```

Regrava `figuras/*.pdf` e `codigo/resultados/metricas.json`. Depois de rodar, confira se os
números das tabelas em `secoes/06-resultados/` continuam batendo com o JSON.

## Como citar e inserir figura

```latex
\cite{butcher2016}                % (Butcher, 2016)
\cite{burden2016,chapra2015}      % (Burden; Faires; Burden, 2016; Chapra; Canale, 2015)
\textcite{butcher2016}            % Butcher (2016)
```

Figura com legenda acima e fonte abaixo:

```latex
\begin{figure}[htbp]
\centering
\caption{Legenda da figura.}
\label{fig:exemplo}
\includegraphics[width=\textwidth]{exemplo.pdf}
\fonte{Elaborado pelo autor (2026).}
\end{figure}
```

Tabelas e algoritmos seguem o mesmo padrão. A figura precisa ser citada no texto
(`Figura~\ref{fig:exemplo}`) antes de aparecer. O que falta confirmar fica marcado com
`\confirmar{...}` (aparece em vermelho no PDF).

## Commits

Conventional Commits em português, com escopo. Exemplo:
`docs(resultados): discute estabilidade com Q = 50`.
