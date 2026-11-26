// =====================================================================
// rlc_experimentos.sci
// Os ensaios do trabalho. Cada função devolve uma struct com os números
// que vão para as tabelas e para o arquivo de resultados.
// =====================================================================

function [nomes, passos, estagios] = metodos()
    // Os três métodos, na ordem das tabelas.
    nomes = ["euler", "heun", "rk4"];
    passos = list(passo_euler, passo_heun, passo_rk4);
    estagios = [1, 2, 4];
endfunction

function v = exp_verificacao(L, C, V0)
    // Confere as fórmulas exatas antes de usá-las como referência.
    v = struct();
    t = linspace(0, 5e-3, 20001);
    nomes = ["sub", "critico", "sobre"];
    Rs = [20, 200, 500];
    for k = 1:3
        p = rlc_parametros(Rs(k), L, C);
        v(nomes(k)) = struct("R", Rs(k), "regime", p.regime, ..
                             "residuo", residuo_degrau(t, V0, p));
    end
    p = rlc_parametros(20, L, C);
    ts = linspace(0, 8e-3, 80001);
    nomes = ["ressonancia", "abaixo", "acima"];
    fatores = [1, 0.7, 1.4];
    for k = 1:3
        w = fatores(k) * p.w0;
        v(nomes(k)) = struct("w", w, "residuo", residuo_seno(ts, V0, w, p));
    end
endfunction

function s = exp_convergencia(p, tf, hs)
    // Erro máximo no degrau para vários passos, ordem observada e ajuste.
    [nomes, passos, estagios] = metodos();
    s = struct();
    for m = 1:3
        es = []; hes = []; evs = []; nfs = [];
        for h = hs
            r = integrar(passos(m), [0; 0], 0, tf, h, p, 1e6);
            e = erro_max(r.y, sol_degrau(r.t, p.V, p));
            hes($ + 1) = r.h;
            es($ + 1) = e(1);
            evs($ + 1) = e(2);
            nfs($ + 1) = r.nf;
        end
        s(nomes(m)) = struct("h", hes', "max_i", es', "max_vc", evs', "avaliacoes", nfs', ..
                             "ordens", ordem_observada(hes, es), ..
                             "ajuste", ajuste_ordem(hes, es));
    end
endfunction

function s = exp_regimes(L, C, V0, tf, h)
    // Erro máximo da corrente com o mesmo passo nos três regimes.
    [nomes, passos, estagios] = metodos();
    s = struct();
    regimes = ["subamortecido", "critico", "sobreamortecido"];
    Rs = [20, 200, 500];
    for k = 1:3
        p = rlc_parametros(Rs(k), L, C);
        p.V = V0;
        item = struct("R", Rs(k));
        for m = 1:3
            r = integrar(passos(m), [0; 0], 0, tf, h, p, 1e6);
            e = erro_max(r.y, sol_degrau(r.t, V0, p));
            item(nomes(m)) = struct("max_i", e(1), "max_vc", e(2));
        end
        s(regimes(k)) = item;
    end
endfunction

function s = exp_estabilidade(L, C, V0, tf, hs)
    // Limite teórico de estabilidade e pico de |i| para passos crescentes,
    // nos circuitos com Q = 5 (R = 20) e Q = 50 (R = 2).
    [nomes, passos, estagios] = metodos();
    s = struct();
    casos = ["Q5", "Q50"];
    Rs = [20, 2];
    for k = 1:2
        p = rlc_parametros(Rs(k), L, C);
        p.V = V0;
        lim = struct();
        item = struct("R", Rs(k), "Q", p.Q);
        for m = 1:3
            lim(nomes(m)) = limite_estab(nomes(m), p);
            hes = []; picos = []; erros = [];
            for h = hs
                r = integrar(passos(m), [0; 0], 0, tf, h, p, 1e3);
                hes($ + 1) = r.h;
                if r.estavel then
                    picos($ + 1) = max(abs(r.y(1, :)));
                    e = erro_max(r.y, sol_degrau(r.t, V0, p));
                    erros($ + 1) = e(1);
                else
                    picos($ + 1) = %nan;
                    erros($ + 1) = %nan;
                end
            end
            item(nomes(m)) = struct("h", hes', "pico_i", picos', "max_i", erros');
        end
        item.h_estabilidade = lim;
        s(casos(k)) = item;
    end
endfunction

function s = exp_senoidal(p, tf, h)
    // Senoide na ressonância e fora dela: fasor e erro de cada método.
    [nomes, passos, estagios] = metodos();
    s = struct();
    casos = ["ressonancia", "abaixo", "acima"];
    fatores = [1, 0.7, 1.4];
    for k = 1:3
        q = p;
        q.senoidal = %t;
        q.w = fatores(k) * p.w0;
        [amp, fase] = fasor_corrente(q.V, q.w, q);
        item = struct("w", q.w, "fator", fatores(k), "amplitude", amp, "fase", fase);
        for m = 1:3
            r = integrar(passos(m), [0; 0], 0, tf, h, q, 1e6);
            e = erro_max(r.y, sol_seno(r.t, q.V, q.w, q));
            item(nomes(m)) = struct("max_i", e(1), "max_vc", e(2));
        end
        s(casos(k)) = item;
    end
endfunction

function s = exp_custo(conv, tf, alvos)
    // Avaliações de f estimadas para cada erro alvo.
    [nomes, passos, estagios] = metodos();
    s = struct();
    for m = 1:3
        s(nomes(m)) = custo_para_erro(conv(nomes(m)).ajuste, estagios(m), tf, alvos);
    end
endfunction

// ---------------------------------------------------------------------
// Comparação extra com o integrador pronto do Scilab (comando ode).
// Os métodos do trabalho continuam escritos à mão; o ode entra só como
// outra técnica numérica para comparar erro e custo (Aquino, 2020, 10.7).
// ---------------------------------------------------------------------

function dy = campo_contado(t, y, p)
    // Mesmo campo do circuito, contando quantas vezes o ode o chama.
    global NF_ODE
    NF_ODE = NF_ODE + 1;
    dy = rlc_campo(t, y, p);
endfunction

function s = exp_ode(p, tf, h, tolerancias)
    // Degrau subamortecido nos mesmos instantes do RK4 com passo h.
    global NF_ODE
    n = round(tf / h);
    t = (tf / n) * (0:n);
    t(n + 1) = tf;
    ref = sol_degrau(t, p.V, p);
    s = struct();
    tipos = ["rk", "lsoda"];
    for k = 1:2
        erros = []; nfs = [];
        for rtol = tolerancias
            NF_ODE = 0;
            atol = 1e-3 * rtol;
            if tipos(k) == "rk" then
                y = ode("rk", [0; 0], 0, t, rtol, atol, list(campo_contado, p));
            else
                y = ode([0; 0], 0, t, rtol, atol, list(campo_contado, p));
            end
            e = erro_max([[0; 0], y(:, 2:$)], ref);
            erros($ + 1) = e(1);
            nfs($ + 1) = NF_ODE;
        end
        s(tipos(k)) = struct("rtol", tolerancias, "max_i", erros', "avaliacoes", nfs');
    end
    r = integrar(passo_rk4, [0; 0], 0, tf, h, p, 1e6);
    e = erro_max(r.y, ref);
    s.rk4_passo_fixo = struct("h", r.h, "max_i", e(1), "avaliacoes", r.nf);
endfunction
