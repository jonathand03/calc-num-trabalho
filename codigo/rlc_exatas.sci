// =====================================================================
// rlc_exatas.sci
// Soluções exatas (analíticas) do RLC série com i(0) = vC(0) = 0.
// Servem de referência para medir o erro dos métodos numéricos.
// Cada função recebe um vetor linha de tempos t e devolve uma matriz
// 2 x N: linha 1 = corrente, linha 2 = tensão no capacitor.
// =====================================================================

function y = sol_degrau(t, V0, p)
    // Resposta a um degrau de V0 volts nos três regimes.
    a = p.alpha;
    e = exp(-a * t);
    L = p.L;
    select p.regime
    case "subamortecido" then
        wd = p.wd;
        i = (V0 / (wd * L)) * e .* sin(wd * t);
        vc = V0 * (1 - e .* (cos(wd * t) + (a / wd) * sin(wd * t)));
    case "sobreamortecido" then
        b = p.beta;
        i = (V0 / (b * L)) * e .* sinh(b * t);
        vc = V0 * (1 - e .* (cosh(b * t) + (a / b) * sinh(b * t)));
    else
        i = (V0 / L) * t .* e;
        vc = V0 * (1 - e .* (1 + a * t));
    end
    y = [i; vc];
endfunction

function di = derivada_degrau(t, V0, p)
    // di/dt exata do degrau. Usada só para conferir a fórmula.
    a = p.alpha;
    e = exp(-a * t);
    select p.regime
    case "subamortecido" then
        wd = p.wd;
        di = (V0 / (wd * p.L)) * e .* (-a * sin(wd * t) + wd * cos(wd * t));
    case "sobreamortecido" then
        b = p.beta;
        di = (V0 / (b * p.L)) * e .* (-a * sinh(b * t) + b * cosh(b * t));
    else
        di = (V0 / p.L) * e .* (1 - a * t);
    end
endfunction

function r = residuo_degrau(t, V0, p)
    // Confere a fórmula pela lei das malhas: V0 - (R i + L di/dt + vC)
    // deve ser zero. Devolve o maior desvio em volts.
    y = sol_degrau(t, V0, p);
    di = derivada_degrau(t, V0, p);
    r = max(abs(V0 - (p.R * y(1, :) + p.L * di + y(2, :))));
endfunction

function [ressonante, D, E, A, B] = coef_seno(Vm, w, p)
    // Coeficientes da corrente para vs = Vm sen(w t), regime subamortecido.
    // Parte permanente: D cos(w t) + E sen(w t), obtida de um sistema 2x2.
    // Parte transitória: exp(-alpha t) (A cos(wd t) + B sen(wd t)),
    // com A e B fixados pelas condições iniciais nulas.
    G = 1 / p.C - p.L * w * w;
    ressonante = abs(G) <= 1e-9 * p.w0^2 * p.L;
    D = 0; E = 0; A = 0; B = 0;
    if ~ressonante then
        M = [G, p.R * w; -p.R * w, G];
        x = M \ [Vm * w; 0];
        D = x(1);
        E = x(2);
        A = -D;
        B = -(E * w + p.alpha * D) / p.wd;
    end
endfunction

function y = sol_seno(t, Vm, w, p)
    // Resposta completa a vs = Vm sen(w t). A tensão no capacitor sai da
    // lei das malhas: vC = vs - R i - L di/dt.
    a = p.alpha;
    wd = p.wd;
    R = p.R;
    [ressonante, D, E, A, B] = coef_seno(Vm, w, p);
    if ressonante then
        // Na ressonância, Z = R e a fórmula fica mais simples.
        i = (Vm / R) * (sin(w * t) - (w / wd) * exp(-a * t) .* sin(wd * t));
        di = (Vm / R) * w * (cos(w * t) - exp(-a * t) .* cos(wd * t));
        di = di + (Vm / R) * (w * a / wd) * exp(-a * t) .* sin(wd * t);
    else
        ex = exp(-a * t);
        i = D * cos(w * t) + E * sin(w * t) + ex .* (A * cos(wd * t) + B * sin(wd * t));
        di_p = -D * w * sin(w * t) + E * w * cos(w * t);
        di_h = -a * ex .* (A * cos(wd * t) + B * sin(wd * t));
        di_h = di_h + ex .* (-A * wd * sin(wd * t) + B * wd * cos(wd * t));
        di = di_p + di_h;
    end
    vc = Vm * sin(w * t) - R * i - p.L * di;
    y = [i; vc];
endfunction

function r = residuo_seno(t, Vm, w, p)
    // Confere a fórmula na EDO da corrente:
    // L i'' + R i' + i/C deve ser igual a d(vs)/dt = Vm w cos(w t).
    a = p.alpha;
    wd = p.wd;
    [ressonante, D, E, A, B] = coef_seno(Vm, w, p);
    e = exp(-a * t);
    sh = sin(wd * t);
    ch = cos(wd * t);
    if ressonante then
        k = Vm / p.R;
        s = sin(w * t);
        c = cos(w * t);
        i = k * (s - (w / wd) * e .* sh);
        di = k * (w * c - (w / wd) * e .* (-a * sh + wd * ch));
        Tp = e .* ((a * a - wd * wd) * sh - 2 * a * wd * ch);
        ddi = k * (-(w^2) * s - (w / wd) * Tp);
    else
        cw = cos(w * t);
        sw = sin(w * t);
        i = D * cw + E * sw + e .* (A * ch + B * sh);
        uh = A * ch + B * sh;
        duh = -A * wd * sh + B * wd * ch;
        dduh = -A * wd * wd * ch - B * wd * wd * sh;
        di = -D * w * sw + E * w * cw + (-a * e .* uh + e .* duh);
        ddi = -D * w * w * cw - E * w * w * sw ..
              + (a * a * e .* uh - 2 * a * e .* duh + e .* dduh);
    end
    r = max(abs(p.L * ddi + p.R * di + i / p.C - Vm * w * cos(w * t)));
endfunction

function [amp, fase] = fasor_corrente(Vm, w, p)
    // Amplitude e fase da corrente em regime permanente senoidal:
    // I = Vm / Z, com Z = R + j (w L - 1/(w C)).
    Z = p.R + %i * (w * p.L - 1 / (w * p.C));
    I = Vm / Z;
    amp = abs(I);
    fase = atan(imag(I), real(I));
endfunction
