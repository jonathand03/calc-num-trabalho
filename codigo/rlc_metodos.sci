// =====================================================================
// rlc_metodos.sci
// Modelo do circuito RLC série e os três métodos de passo único:
// Euler, Euler melhorado (Heun) e Runge-Kutta clássico de quarta ordem.
// As fórmulas dos passos seguem Aquino (2020, seções 10.1, 10.2 e 10.3.2).
// O estado é o vetor coluna y = [i; vC]: corrente e tensão no capacitor.
// =====================================================================

function p = rlc_parametros(R, L, C)
    // Grandezas derivadas do circuito: amortecimento (alpha), pulsação
    // natural (w0), pulsação amortecida (wd) e fator de qualidade (Q).
    p = struct();
    p.R = R;
    p.L = L;
    p.C = C;
    p.w0 = 1 / sqrt(L * C);
    p.f0 = p.w0 / (2 * %pi);
    p.alpha = R / (2 * L);
    p.wd = %nan;
    p.beta = %nan;
    disc = p.alpha^2 - p.w0^2;
    if disc < -1e-18 then
        p.regime = "subamortecido";
        p.wd = sqrt(p.w0^2 - p.alpha^2);
    elseif disc > 1e-18 then
        p.regime = "sobreamortecido";
        p.beta = sqrt(disc);
    else
        p.regime = "critico";
    end
    p.zeta = p.alpha / p.w0;
    p.Q = p.w0 * L / R;
    // Raízes de lambda^2 + 2 alpha lambda + w0^2 = 0 (autovalores do sistema).
    raiz = sqrt(p.alpha^2 - p.w0^2);
    p.lambda = [-p.alpha + raiz; -p.alpha - raiz];
    // Fonte: degrau por padrão. A senoide liga "senoidal" e define "w".
    p.V = 5;
    p.senoidal = %f;
    p.w = 0;
endfunction

function v = fonte_vs(t, p)
    // Tensão da fonte no instante t: degrau de valor V ou V sen(w t).
    if p.senoidal then
        v = p.V * sin(p.w * t);
    else
        v = p.V;
    end
endfunction

function dy = rlc_campo(t, y, p)
    // Lado direito de y' = f(t, y).
    // Lei das malhas: L di/dt = vs - R i - vC.  Capacitor: C dvC/dt = i.
    di = (fonte_vs(t, p) - p.R * y(1) - y(2)) / p.L;
    dvc = y(1) / p.C;
    dy = [di; dvc];
endfunction

function [y1, nf] = passo_euler(t, y, h, p)
    // Método de Euler (Aquino, 2020, eq. 10.4):
    //     y(k+1) = y(k) + h f(t(k), y(k))
    y1 = y + h * rlc_campo(t, y, p);
    nf = 1;
endfunction

function [y1, nf] = passo_heun(t, y, h, p)
    // Euler melhorado ou Heun (Aquino, 2020, eq. 10.6):
    //     K1 = y(k) + h f(t(k), y(k))                          (previsão)
    //     y(k+1) = y(k) + h/2 [f(t(k), y(k)) + f(t(k)+h, K1)]  (correção)
    // f(t(k), y(k)) é calculada uma vez só e reaproveitada.
    f1 = rlc_campo(t, y, p);
    K1 = y + h * f1;
    y1 = y + h * (f1 + rlc_campo(t + h, K1, p)) / 2;
    nf = 2;
endfunction

function [y1, nf] = passo_rk4(t, y, h, p)
    // Runge-Kutta de quarta ordem (Aquino, 2020, seção 10.3.2):
    //     K1 = h f(t(k), y(k))
    //     K2 = h f(t(k) + h/2, y(k) + K1/2)
    //     K3 = h f(t(k) + h/2, y(k) + K2/2)
    //     K4 = h f(t(k) + h, y(k) + K3)
    //     y(k+1) = y(k) + (K1 + 2 K2 + 2 K3 + K4)/6
    K1 = h * rlc_campo(t, y, p);
    K2 = h * rlc_campo(t + h/2, y + K1/2, p);
    K3 = h * rlc_campo(t + h/2, y + K2/2, p);
    K4 = h * rlc_campo(t + h, y + K3, p);
    y1 = y + (K1 + 2*K2 + 2*K3 + K4) / 6;
    nf = 4;
endfunction

function r = integrar(passo, y0, t0, tf, h, p, limite)
    // Integra de t0 até tf com passo fixo.
    // O passo é ajustado para que a malha termine exatamente em tf.
    // Critério de parada: chegar a tf, ou o estado deixar de ser finito,
    // ou algum componente passar de "limite" (divergência numérica).
    n = max(round((tf - t0) / h), 1);
    he = (tf - t0) / n;
    t = t0 + he * (0:n);
    t(n + 1) = tf;
    y = %nan * ones(2, n + 1);
    y(:, 1) = y0;
    nf = 0;
    estavel = %t;
    for k = 1:n
        [yn, m] = passo(t(k), y(:, k), he, p);
        nf = nf + m;
        if or(isnan(yn)) | or(isinf(yn)) | max(abs(yn)) > limite then
            estavel = %f;
            break
        end
        y(:, k + 1) = yn;
    end
    r = struct("t", t, "y", y, "nf", nf, "estavel", estavel, ..
               "h", he, "passos", n);
endfunction
