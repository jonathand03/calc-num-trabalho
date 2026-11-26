// =====================================================================
// main.sce
// Simulação do circuito RLC série pelos métodos de Euler, Euler melhorado
// e Runge-Kutta de quarta ordem. Uma execução gera todos os números das
// tabelas (codigo/resultados/metricas_scilab.json) e as figuras (figuras/).
//
// Como rodar, no console do Scilab:
//     exec("codigo/main.sce", -1)
// ou no terminal, na raiz do repositório:
//     scilab-cli -nb -f codigo/main.sce        (sem figuras)
//     scilab -nw -nb -f codigo/main.sce        (com figuras)
// =====================================================================

clear;
funcprot(0);
base = get_absolute_file_path("main.sce");
exec(base + "rlc_metodos.sci", -1);
exec(base + "rlc_exatas.sci", -1);
exec(base + "rlc_analise.sci", -1);
exec(base + "rlc_experimentos.sci", -1);
exec(base + "rlc_figuras.sci", -1);
pasta_fig = base + "../figuras/";
pasta_res = base + "resultados/";
mkdir(pasta_res);

// ---------------------------------------------------------------------
// Entrada de dados
// ---------------------------------------------------------------------
L = 10e-3;          // indutância (H)
C = 1e-6;           // capacitância (F)
V0 = 5;             // amplitude da fonte (V)
R = 20;             // resistência do circuito principal (ohm), Q = 5
tf_conv = 5e-3;     // intervalo dos ensaios de convergência (s)
tf_estab = 6e-3;    // intervalo do ensaio de estabilidade (s)
tf_seno = 8e-3;     // intervalo do ensaio senoidal (s)
hs_conv = [8e-6, 4e-6, 2e-6, 1e-6, 5e-7, 2.5e-7];
hs_estab = [5, 10, 15, 20, 25, 40, 80, 120, 200, 280] * 1e-6;
alvos = [1e-3, 1e-4, 1e-6, 1e-8];

p = rlc_parametros(R, L, C);
p.V = V0;

// ---------------------------------------------------------------------
// Ensaios
// ---------------------------------------------------------------------
mprintf("Conferindo as soluções exatas...\n");
verif = exp_verificacao(L, C, V0);
for nome = fieldnames(verif)'
    mprintf("  %-12s resíduo = %.3e\n", nome, verif(nome).residuo);
    if verif(nome).residuo > 5e-3 then
        error("Solução exata inconsistente em " + nome);
    end
end

mprintf("Convergência...\n");
conv = exp_convergencia(p, tf_conv, hs_conv);
mprintf("Regimes...\n");
regimes = exp_regimes(L, C, V0, tf_conv, 2e-6);
mprintf("Estabilidade...\n");
estab = exp_estabilidade(L, C, V0, tf_estab, hs_estab);
mprintf("Senoide...\n");
seno = exp_senoidal(p, tf_seno, 1e-6);
mprintf("Custo...\n");
custo = exp_custo(conv, tf_conv, alvos);
mprintf("Comparação com o ode do Scilab...\n");
cmp_ode = exp_ode(p, tf_conv, 2e-6, [1e-6, 1e-10]);

// ---------------------------------------------------------------------
// Resultados
// ---------------------------------------------------------------------
circuito = struct();
nomes_c = ["principal", "Q50", "critico", "sobre"];
Rs_c = [20, 2, 200, 500];
for k = 1:4
    q = rlc_parametros(Rs_c(k), L, C);
    circuito(nomes_c(k)) = struct("R", q.R, "w0", q.w0, "f0", q.f0, "alpha", q.alpha, ..
        "wd", q.wd, "beta", q.beta, "zeta", q.zeta, "Q", q.Q, "regime", q.regime, ..
        "periodo_amortecido", 2 * %pi / q.wd, "lambda_real", real(q.lambda)', ..
        "lambda_imag", imag(q.lambda)', "RC", q.R * q.C);
end

resultado = struct("base", struct("L", L, "C", C, "V0", V0, "eps_dupla", %eps), ..
    "circuitos", circuito, "verificacao", verif, "convergencia", conv, ..
    "regimes", regimes, "estabilidade", estab, "senoidal", seno, ..
    "custo", custo, "ode", cmp_ode, "versao_scilab", getversion());
gravar_json(resultado, pasta_res + "metricas_scilab.json");
mprintf("Resultados em %s\n", pasta_res + "metricas_scilab.json");

// ---------------------------------------------------------------------
// Figuras (precisam do Scilab com interface gráfica, não do scilab-cli)
// ---------------------------------------------------------------------
if getscilabmode() <> "NWNI" then
    mprintf("Gerando figuras...\n");
    fig_convergencia(conv, pasta_fig + "convergencia.pdf");
    fig_degrau(p, tf_estab, 15e-6, pasta_fig + "resposta_degrau.pdf");
    fig_erro_temporal(p, tf_conv, [8e-6, 2e-6], pasta_fig + "erro_temporal.pdf");
    fig_estabilidade(estab, pasta_fig + "estabilidade.pdf");
    fig_regimes(L, C, V0, tf_conv, 2e-6, pasta_fig + "regimes.pdf");
    fig_senoidal(p, tf_seno, 1e-6, pasta_fig + "senoidal.pdf");
else
    mprintf("Modo sem gráficos: figuras não geradas.\n");
end

// ---------------------------------------------------------------------
// Resumo no console
// ---------------------------------------------------------------------
mprintf("\nPeríodo amortecido (R = 20): %.4f us\n", 2 * %pi / p.wd * 1e6);
mprintf("\nLimites de estabilidade (us):\n");
for c = ["Q5", "Q50"]
    lim = estab(c).h_estabilidade;
    mprintf("  %-4s Euler %.2f  Heun %.2f  RK4 %.2f\n", c, ..
            lim.euler * 1e6, lim.heun * 1e6, lim.rk4 * 1e6);
end
mprintf("\nAjuste de ordem p:\n");
for m = ["euler", "heun", "rk4"]
    mprintf("  %-6s p = %.4f\n", m, conv(m).ajuste.p);
end
mprintf("\nComparação com o ode (degrau, h de saída = 2 us):\n");
for tp = ["rk", "lsoda"]
    o = cmp_ode(tp);
    for k = 1:size(o.rtol, "*")
        mprintf("  %-6s rtol %.0e  erro %.3e A  avaliações %d\n", ..
                tp, o.rtol(k), o.max_i(k), o.avaliacoes(k));
    end
end
mprintf("  RK4 próprio  erro %.3e A  avaliações %d\n", ..
        cmp_ode.rk4_passo_fixo.max_i, cmp_ode.rk4_passo_fixo.avaliacoes);
