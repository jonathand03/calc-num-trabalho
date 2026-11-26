// =====================================================================
// rlc_figuras.sci
// Gráficos do trabalho, exportados em PDF para a pasta figuras/.
// =====================================================================

function c = cor_metodo(nome)
    // Cores fixas por curva, legíveis em impressão e para daltonismo.
    select nome
    case "analitica" then
        c = color(26, 26, 26);
    case "euler" then
        c = color(192, 57, 43);
    case "heun" then
        c = color(36, 113, 163);
    case "rk4" then
        c = color(30, 132, 73);
    else
        c = color(127, 140, 141);
    end
endfunction

function r = rotulo(nome)
    select nome
    case "euler" then
        r = "Euler";
    case "heun" then
        r = "Euler melhorado";
    case "rk4" then
        r = "RK4";
    else
        r = "Solução analítica";
    end
endfunction

function l = linha(x, y, nome, espessura, estilo)
    // Desenha uma curva e devolve o objeto, para montar a legenda.
    plot(x, y);
    e = gce();
    l = e.children(1);
    l.foreground = cor_metodo(nome);
    l.thickness = espessura;
    l.line_style = estilo;
endfunction

function f = nova_figura(largura, altura)
    f = scf();
    clf(f);
    f.figure_size = [largura, altura];
    f.background = color(255, 255, 255);
endfunction

function eixos(titulo, rx, ry)
    // Título, rótulos e fonte dos eixos atuais.
    a = gca();
    a.font_size = 2;
    a.title.text = titulo;
    a.title.font_size = 3;
    a.x_label.text = rx;
    a.x_label.font_size = 2;
    a.y_label.text = ry;
    a.y_label.font_size = 2;
    a.grid = [color(220, 220, 220), color(220, 220, 220)];
    a.background = color(255, 255, 255);
endfunction

function salvar(f, arquivo)
    xs2pdf(f, arquivo);
    delete(f);
endfunction

function fig_degrau(p, tf, h, arquivo)
    // Corrente e tensão no capacitor para o degrau, três métodos e exata.
    [nomes, passos, estagios] = metodos();
    tfino = linspace(0, tf, 4000);
    ref = sol_degrau(tfino, p.V, p);
    for m = 1:3
        r(m) = integrar(passos(m), [0; 0], 0, tf, h, p, 1e6);
    end
    f = nova_figura(740, 620);
    escalas = [1e3, 1];
    titulos = ["Resposta ao degrau de 5 V, passo h = " ..
               + msprintf("%.0f", h * 1e6) + " µs", ""];
    ylab = ["Corrente (mA)", "Tensão no capacitor (V)"];
    for k = 1:2
        subplot(2, 1, k);
        hs = linha(tfino * 1e3, ref(k, :) * escalas(k), "analitica", 3, 1);
        for m = 1:3
            hs($ + 1) = linha(r(m).t * 1e3, r(m).y(k, :) * escalas(k), nomes(m), 1.5, 1);
        end
        eixos(titulos(k), "Tempo (ms)", ylab(k));
        if k == 1 then
            legend(hs, ["Solução analítica", "Euler", "Euler melhorado", "RK4"], ..
                   "in_upper_right");
        end
    end
    salvar(f, arquivo);
endfunction

function fig_erro_temporal(p, tf, hs, arquivo)
    // Erro absoluto da corrente ao longo do tempo, para dois passos.
    [nomes, passos, estagios] = metodos();
    f = nova_figura(760, 360);
    faixas = [1e-14, 1e-1; 1e-17, 1e-2];
    for k = 1:2
        subplot(1, 2, k);
        hl = [];
        for m = 1:3
            r = integrar(passos(m), [0; 0], 0, tf, hs(k), p, 1e6);
            ref = sol_degrau(r.t, p.V, p);
            erro = max(abs(r.y(1, 2:$) - ref(1, 2:$)), 1e-18);
            hl($ + 1) = linha(r.t(2:$) * 1e3, erro, nomes(m), 1.5, 1);
        end
        a = gca();
        a.log_flags = "nl";
        // Faixa vertical com folga embaixo, onde fica a legenda.
        a.data_bounds = [0, faixas(k, 1); tf * 1e3, faixas(k, 2)];
        a.tight_limits = "on";
        eixos("h = " + strsubst(msprintf("%.1f", hs(k) * 1e6), ".", ",") + " µs", ..
              "Tempo (ms)", "Erro absoluto em i (A)");
        legend(hl, ["Euler", "Euler melhorado", "RK4"], "in_lower_right");
    end
    salvar(f, arquivo);
endfunction

function fig_convergencia(conv, arquivo)
    // Erro máximo contra o passo, em escala log-log, com retas de
    // inclinação 1, 2 e 4 ancoradas no maior passo de cada método.
    [nomes, passos, estagios] = metodos();
    marcas = [9, 11, 5];
    ordens = [1, 2, 4];
    estilos = [3, 2, 4];
    f = nova_figura(720, 440);
    hl = [];
    href = [2.5e-7, 8e-6];
    for m = 1:3
        c = conv(nomes(m));
        l = linha(c.h, c.max_i, nomes(m), 2, 1);
        l.mark_mode = "on";
        l.mark_style = marcas(m);
        l.mark_size = 6;
        l.mark_foreground = cor_metodo(nomes(m));
        l.mark_background = cor_metodo(nomes(m));
        hl($ + 1) = l;
        linha(href, c.max_i(1) * (href / c.h(1)) .^ ordens(m), nomes(m), 1, estilos(m));
    end
    a = gca();
        a.log_flags = "ll";
    eixos("Convergência no degrau subamortecido", "Passo h (s)", "Erro máximo em i (A)");
    legend(hl, ["Euler", "Euler melhorado", "RK4"], "in_lower_right");
    salvar(f, arquivo);
endfunction

function fig_estabilidade(estab, arquivo)
    // Pico de |i| contra o passo. As retas verticais pontilhadas são os
    // limites de estabilidade de cada método.
    [nomes, passos, estagios] = metodos();
    casos = ["Q5", "Q50"];
    f = nova_figura(760, 380);
    for k = 1:2
        subplot(1, 2, k);
        b = estab(casos(k));
        hl = [];
        for m = 1:3
            s = b(nomes(m));
            picos = min(s.pico_i, 50);
            l = linha(s.h * 1e6, picos, nomes(m), 1.5, 1);
            l.mark_mode = "on";
            l.mark_style = 9;
            l.mark_size = 4;
            l.mark_foreground = cor_metodo(nomes(m));
            hl($ + 1) = l;
        end
        for m = 1:3
            hc = b.h_estabilidade(nomes(m)) * 1e6;
            linha([hc, hc], [1e-2, 50], nomes(m), 1, 3);
        end
        a = gca();
        a.log_flags = "nl";
        // Folga acima de 50 A (teto do gráfico) para a legenda.
        a.data_bounds = [0, 1e-2; 300, 1e5];
        a.tight_limits = "on";
        titulo = "R = " + msprintf("%.0f", b.R) + " Ω, Q = " + msprintf("%.0f", b.Q);
        eixos(titulo, "Passo h (µs)", "Pico de |i| no intervalo (A)");
        legend(hl, ["Euler", "Euler melhorado", "RK4"], "in_upper_right");
    end
    salvar(f, arquivo);
endfunction

function fig_regimes(L, C, V0, tf, h, arquivo)
    // Corrente do degrau nos três regimes: exata e RK4.
    Rs = [20, 200, 500];
    titulos = ["Subamortecido, R = 20 Ω", "Crítico, R = 200 Ω", ..
               "Sobreamortecido, R = 500 Ω"];
    f = nova_figura(740, 740);
    for k = 1:3
        subplot(3, 1, k);
        p = rlc_parametros(Rs(k), L, C);
        p.V = V0;
        t = linspace(0, tf, 3000);
        ref = sol_degrau(t, V0, p);
        r = integrar(passo_rk4, [0; 0], 0, tf, h, p, 1e6);
        h1 = linha(t * 1e3, ref(1, :) * 1e3, "analitica", 3, 1);
        h2 = linha(r.t * 1e3, r.y(1, :) * 1e3, "rk4", 1.5, 1);
        rx = "";
        if k == 3 then
            rx = "Tempo (ms)";
        end
        eixos(titulos(k), rx, "i (mA)");
        legend([h1, h2], ["Analítica", "RK4"], "in_upper_right");
    end
    salvar(f, arquivo);
endfunction

function fig_senoidal(p, tf, h, arquivo)
    // Corrente com fonte senoidal, na ressonância e em 0,7 w0.
    fatores = [1, 0.7];
    titulos = ["Ressonância, ω = ω₀", "Fora da ressonância, ω = 0,7 ω₀"];
    ylim = [300, 100];
    f = nova_figura(740, 600);
    for k = 1:2
        subplot(2, 1, k);
        q = p;
        q.senoidal = %t;
        q.w = fatores(k) * p.w0;
        t = linspace(0, tf, 5000);
        ref = sol_seno(t, q.V, q.w, q);
        r = integrar(passo_rk4, [0; 0], 0, tf, h, q, 1e6);
        [amp, fase] = fasor_corrente(q.V, q.w, q);
        h1 = linha(t * 1e3, ref(1, :) * 1e3, "analitica", 3, 1);
        h2 = linha(r.t * 1e3, r.y(1, :) * 1e3, "rk4", 1.5, 1);
        h3 = linha([0, tf * 1e3], [amp, amp] * 1e3, "cinza", 1, 2);
        linha([0, tf * 1e3], -[amp, amp] * 1e3, "cinza", 1, 2);
        rx = "";
        if k == 2 then
            rx = "Tempo (ms)";
        end
        // Folga acima das curvas para a legenda.
        a = gca();
        a.data_bounds = [0, -ylim(k); tf * 1e3, 1.8 * ylim(k)];
        a.tight_limits = "on";
        eixos(titulos(k), rx, "Corrente (mA)");
        legend([h1, h2, h3], ["Analítica", "RK4", "Amplitude permanente"], ..
               "in_upper_right");
    end
    salvar(f, arquivo);
endfunction
