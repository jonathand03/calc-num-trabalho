#!/usr/bin/env python3
"""
Simulação de um circuito RLC série pelos métodos de Euler, Euler melhorado
(Heun) e Runge-Kutta clássico de quarta ordem.

Os métodos são implementados neste arquivo. A solução analítica do circuito
linear é usada como referência. NumPy entra apenas como estrutura de vetores;
Matplotlib gera as figuras.

Executar a partir de qualquer diretório:
    python3 codigo/rlc_runge_kutta.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figuras"
RES = Path(__file__).resolve().parent / "resultados"

# Paleta estável para daltonismo e impressão.
COR = {
    "analitica": "#1a1a1a",
    "euler": "#c0392b",
    "heun": "#2471a3",
    "rk4": "#1e8449",
    "pulso": "#6c3483",
    "float32": "#b9770e",
}


def configurar_estilo() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.35,
            "grid.linestyle": "--",
            "figure.dpi": 140,
            "savefig.dpi": 200,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.axisbelow": True,
        }
    )


# ---------------------------------------------------------------------------
# Modelo do circuito
# ---------------------------------------------------------------------------

def parametros_derivados(R: float, L: float, C: float) -> dict:
    """Grandezas do RLC série: amortecimento, frequência natural e fator Q."""
    w0 = 1.0 / np.sqrt(L * C)
    alpha = R / (2.0 * L)
    discriminante = alpha**2 - w0**2
    if discriminante < -1e-18:
        regime = "subamortecido"
        wd = float(np.sqrt(w0**2 - alpha**2))
        beta = None
    elif discriminante > 1e-18:
        regime = "sobreamortecido"
        wd = None
        beta = float(np.sqrt(discriminante))
    else:
        regime = "critico"
        wd = None
        beta = None
    Q = w0 * L / R
    autovalores = autovalores_sistema(R, L, C)
    return {
        "R": R,
        "L": L,
        "C": C,
        "w0": float(w0),
        "f0": float(w0 / (2.0 * np.pi)),
        "alpha": float(alpha),
        "wd": None if wd is None else float(wd),
        "beta": beta,
        "zeta": float(alpha / w0),
        "Q": float(Q),
        "regime": regime,
        "autovalores_real": [float(np.real(l)) for l in autovalores],
        "autovalores_imag": [float(np.imag(l)) for l in autovalores],
        "h_crit_euler_formula": float(R * C) if regime == "subamortecido" else None,
    }


def autovalores_sistema(R: float, L: float, C: float) -> np.ndarray:
    """Autovalores de y' = A y para y = [i, v_C]."""
    w0_2 = 1.0 / (L * C)
    alpha = R / (2.0 * L)
    # lambda^2 + 2 alpha lambda + omega0^2 = 0
    return np.roots([1.0, 2.0 * alpha, w0_2]).astype(complex)


def fazer_campo(R: float, L: float, C: float, fonte, dtype=np.float64):
    """Campo f(t, y) do sistema de estados, no tipo numérico pedido."""
    Rd, Ld, Cd = dtype(R), dtype(L), dtype(C)

    def f(t, y):
        tt = dtype(t) if not isinstance(t, np.floating) else t.astype(dtype, copy=False)
        i = y[0]
        vc = y[1]
        di = (dtype(fonte(float(tt))) - Rd * i - vc) / Ld
        dvc = i / Cd
        return np.array([di, dvc], dtype=dtype)

    return f


def fonte_degrau(V0: float):
    def vs(_t: float) -> float:
        return V0

    return vs


def fonte_seno(Vm: float, w: float):
    def vs(t: float) -> float:
        return Vm * np.sin(w * t)

    return vs


def fonte_pulso(Vm: float, t0: float, sigma: float):
    def vs(t: float) -> float:
        return Vm * np.exp(-0.5 * ((t - t0) / sigma) ** 2)

    return vs


# ---------------------------------------------------------------------------
# Soluções analíticas (condições iniciais nulas, salvo indicação)
# ---------------------------------------------------------------------------

def solucao_degrau(t, V0, R, L, C):
    """Corrente e tensão no capacitor para degrau V0, i(0)=0, v_C(0)=0."""
    info = parametros_derivados(R, L, C)
    a = info["alpha"]
    w0 = info["w0"]
    t = np.asarray(t, dtype=float)
    e = np.exp(-a * t)

    if info["regime"] == "subamortecido":
        wd = info["wd"]
        i = (V0 / (wd * L)) * e * np.sin(wd * t)
        vc = V0 * (1.0 - e * (np.cos(wd * t) + (a / wd) * np.sin(wd * t)))
    elif info["regime"] == "sobreamortecido":
        b = info["beta"]
        i = (V0 / (b * L)) * e * np.sinh(b * t)
        vc = V0 * (1.0 - e * (np.cosh(b * t) + (a / b) * np.sinh(b * t)))
    else:
        i = (V0 / L) * t * e
        vc = V0 * (1.0 - e * (1.0 + a * t))
    return np.column_stack([i, vc])


def solucao_senoidal(t, Vm, w, R, L, C):
    """Resposta completa a v_s = Vm sen(omega t), com i(0) = v_C(0) = 0.

    A tensão no capacitor é recuperada pela lei de Kirchhoff, a partir da
    corrente analítica e da sua derivada.
    """
    info = parametros_derivados(R, L, C)
    if info["regime"] != "subamortecido":
        raise ValueError("A forma fechada senoidal deste script exige regime subamortecido.")
    a = info["alpha"]
    wd = info["wd"]
    w0_2 = info["w0"] ** 2
    t = np.asarray(t, dtype=float)

    G = 1.0 / C - L * w * w
    # Na ressonância o sistema linear para os coeficientes fica singular.
    if abs(G) <= 1e-9 * w0_2 * L:
        i = (Vm / R) * (np.sin(w * t) - (w / wd) * np.exp(-a * t) * np.sin(wd * t))
        di = (Vm / R) * w * (np.cos(w * t) - np.exp(-a * t) * np.cos(wd * t))
        di = di + (Vm / R) * (w * a / wd) * np.exp(-a * t) * np.sin(wd * t)
    else:
        matriz = np.array([[G, R * w], [-R * w, G]], dtype=float)
        D, E = np.linalg.solve(matriz, np.array([Vm * w, 0.0]))
        A = -D
        B = -(E * w + a * D) / wd
        expa = np.exp(-a * t)
        i_p = D * np.cos(w * t) + E * np.sin(w * t)
        i_h = expa * (A * np.cos(wd * t) + B * np.sin(wd * t))
        i = i_p + i_h
        di_p = -D * w * np.sin(w * t) + E * w * np.cos(w * t)
        di_h = -a * expa * (A * np.cos(wd * t) + B * np.sin(wd * t))
        di_h = di_h + expa * (-A * wd * np.sin(wd * t) + B * wd * np.cos(wd * t))
        di = di_p + di_h

    vs = Vm * np.sin(w * t)
    vc = vs - R * i - L * di
    return np.column_stack([i, vc])


def fasor_corrente(Vm, w, R, L, C):
    """Amplitude e fase da corrente em regime permanente senoidal."""
    Z = R + 1j * (w * L - 1.0 / (w * C))
    I = Vm / Z
    return float(np.abs(I)), float(np.angle(I))


def _derivada_degrau(t, V0, R, L, C):
    """di/dt analítica do degrau, usada só na verificação da forma fechada."""
    info = parametros_derivados(R, L, C)
    a = info["alpha"]
    t = np.asarray(t, dtype=float)
    e = np.exp(-a * t)
    if info["regime"] == "subamortecido":
        wd = info["wd"]
        k = V0 / (wd * L)
        return k * e * (-a * np.sin(wd * t) + wd * np.cos(wd * t))
    if info["regime"] == "sobreamortecido":
        b = info["beta"]
        k = V0 / (b * L)
        return k * e * (-a * np.sinh(b * t) + b * np.cosh(b * t))
    return (V0 / L) * e * (1.0 - a * t)


def residuo_degrau_exato(t, V0, R, L, C) -> float:
    """Resíduo de KVL com a derivada analítica da corrente."""
    y = solucao_degrau(t, V0, R, L, C)
    di = _derivada_degrau(t, V0, R, L, C)
    kvl = V0 - (R * y[:, 0] + L * di + y[:, 1])
    return float(np.max(np.abs(kvl)))


def residuo_senoidal(t, Vm, w, R, L, C) -> float:
    """Resíduo da EDO da corrente, com derivadas analíticas da forma fechada.

    L i'' + R i' + i/C deve reproduzir d v_s / dt = Vm omega cos(omega t).
    """
    info = parametros_derivados(R, L, C)
    a = info["alpha"]
    wd = info["wd"]
    t = np.asarray(t, dtype=float)
    G = 1.0 / C - L * w * w
    if abs(G) <= 1e-9 * info["w0"] ** 2 * L:
        # i = (Vm/R) [sen(omega t) - (omega/wd) exp(-a t) sen(wd t)]
        k = Vm / R
        s, c = np.sin(w * t), np.cos(w * t)
        e = np.exp(-a * t)
        sh, ch = np.sin(wd * t), np.cos(wd * t)
        i = k * (s - (w / wd) * e * sh)
        di = k * (w * c - (w / wd) * e * (-a * sh + wd * ch))
        # derivada do termo transitório T = e ( -a sen + wd cos )
        # T' = -a e (-a sen + wd cos) + e (-a wd cos - wd^2 sen)
        Tp = e * ((a * a - wd * wd) * sh - 2.0 * a * wd * ch)
        ddi = k * (-(w**2) * s - (w / wd) * Tp)
    else:
        matriz = np.array([[G, R * w], [-R * w, G]], dtype=float)
        D, E = np.linalg.solve(matriz, np.array([Vm * w, 0.0]))
        A = -D
        B = -(E * w + a * D) / wd
        e = np.exp(-a * t)
        cw, sw = np.cos(w * t), np.sin(w * t)
        ch, sh = np.cos(wd * t), np.sin(wd * t)
        i = D * cw + E * sw + e * (A * ch + B * sh)
        di_p = -D * w * sw + E * w * cw
        ddi_p = -D * w * w * cw - E * w * w * sw
        uh = A * ch + B * sh
        duh = -A * wd * sh + B * wd * ch
        dduh = -A * wd * wd * ch - B * wd * wd * sh
        di = di_p + (-a * e * uh + e * duh)
        ddi = ddi_p + (a * a * e * uh - 2.0 * a * e * duh + e * dduh)
    edo = L * ddi + R * di + i / C - Vm * w * np.cos(w * t)
    return float(np.max(np.abs(edo)))


# ---------------------------------------------------------------------------
# Integradores de passo fixo e adaptativo
# ---------------------------------------------------------------------------

def passo_euler(f, t, y, h):
    k1 = f(t, y)
    return y + h * k1, 1


def passo_heun(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h, y + h * k1)
    return y + 0.5 * h * (k1 + k2), 2


def passo_rk4(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4), 4


PASSOS = {
    "euler": passo_euler,
    "heun": passo_heun,
    "rk4": passo_rk4,
}
ESTAGIOS = {"euler": 1, "heun": 2, "rk4": 4}


def malha_tempo(t0: float, tf: float, h: float):
    """Malha uniforme que termina exatamente em tf. O passo é ajustado."""
    n = max(int(round((tf - t0) / h)), 1)
    h_eff = (tf - t0) / n
    t = t0 + h_eff * np.arange(n + 1, dtype=float)
    t[-1] = tf
    return t, h_eff, n


def integrar(passo, f, y0, t0, tf, h, limite=1e6, dtype=np.float64):
    """Integração de passo fixo. Interrompe se a solução divergir."""
    t, h_eff, n = malha_tempo(t0, tf, h)
    y = np.empty((n + 1, len(y0)), dtype=dtype)
    y[0] = np.asarray(y0, dtype=dtype)
    avaliacoes = 0
    estavel = True
    for k in range(n):
        y_novo, nfe = passo(f, t[k], y[k], dtype(h_eff))
        avaliacoes += nfe
        if (not np.isfinite(y_novo).all()) or np.max(np.abs(y_novo)) > limite:
            y[k + 1 :] = np.nan
            estavel = False
            break
        y[k + 1] = y_novo
    return t, y, avaliacoes, estavel, h_eff


def rk4_adaptativo(f, y0, t0, tf, h0, rtol, atol, max_passos=200_000):
    """RK4 com controle de passo por duplicação e extrapolação local.

    O erro local da solução de dois meios passos é estimado por
    ||y_h - y_{h/2}|| / 15, pois o erro local de um método de ordem 4
    escala com h^5 e 2*(h/2)^5 = h^5/16.
    """
    t = float(t0)
    y = np.array(y0, dtype=float)
    h = float(h0)
    tempos = [t]
    estados = [y.copy()]
    passos = []
    avaliacoes = 0
    aceitos = 0
    rejeitados = 0
    seguranca = 0.9
    ordem = 4
    h_min = 1e-12

    while t < tf - 1e-15:
        if aceitos + rejeitados >= max_passos:
            break
        h = min(h, tf - t)
        y_cheio, n1 = passo_rk4(f, t, y, h)
        y_meio, n2 = passo_rk4(f, t, y, 0.5 * h)
        y_metade, n3 = passo_rk4(f, t + 0.5 * h, y_meio, 0.5 * h)
        avaliacoes += n1 + n2 + n3
        delta = y_cheio - y_metade
        escala = atol + rtol * np.maximum(np.abs(y_metade), np.abs(y))
        err = float(np.max(np.abs(delta) / escala) / 15.0)

        if err <= 1.0 or h <= h_min * 1.01:
            y = y_metade
            t = t + h
            tempos.append(t)
            estados.append(y.copy())
            passos.append(h)
            aceitos += 1
            if t >= tf - 1e-15:
                break
            fator = 2.0 if err < 1e-16 else seguranca * err ** (-1.0 / (ordem + 1))
            h = h * min(2.0, max(0.2, fator))
        else:
            rejeitados += 1
            fator = seguranca * err ** (-1.0 / (ordem + 1))
            h = max(h_min, h * min(0.9, max(0.2, fator)))

    return (
        np.array(tempos),
        np.vstack(estados),
        np.array(passos),
        avaliacoes,
        aceitos,
        rejeitados,
    )


def funcao_estabilidade(metodo: str, z: complex) -> complex:
    if metodo == "euler":
        return 1.0 + z
    if metodo == "heun":
        return 1.0 + z + z**2 / 2.0
    if metodo == "rk4":
        return 1.0 + z + z**2 / 2.0 + z**3 / 6.0 + z**4 / 24.0
    raise KeyError(metodo)


def limite_estabilidade(metodo: str, R, L, C, h_teto=5e-3, amostras=40_000) -> float:
    """Maior passo para o qual todos os autovalores caem na região de estabilidade."""
    lambdas = autovalores_sistema(R, L, C)
    hs = np.linspace(h_teto / amostras, h_teto, amostras)
    limite = 0.0
    for h in hs:
        estavel = True
        for lam in lambdas:
            if abs(funcao_estabilidade(metodo, h * lam)) > 1.0 + 1e-9:
                estavel = False
                break
        if estavel:
            limite = float(h)
        else:
            break
    return limite


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------

def erros_contra_referencia(y, y_ref) -> dict:
    mascara = np.isfinite(y).all(axis=1) & np.isfinite(y_ref).all(axis=1)
    if mascara.sum() < 3:
        return {"max_i": None, "rmse_i": None, "max_vc": None, "rmse_vc": None, "max_estado": None}
    di = np.abs(y[mascara, 0] - y_ref[mascara, 0])
    dv = np.abs(y[mascara, 1] - y_ref[mascara, 1])
    return {
        "max_i": float(np.max(di)),
        "rmse_i": float(np.sqrt(np.mean(di**2))),
        "max_vc": float(np.max(dv)),
        "rmse_vc": float(np.sqrt(np.mean(dv**2))),
        "max_estado": float(max(np.max(di), np.max(dv))),
    }


def ordem_observada(hs, erros) -> list:
    pares = []
    for (h1, e1), (h2, e2) in zip(zip(hs, erros), list(zip(hs, erros))[1:]):
        if e1 is None or e2 is None or e1 <= 0 or e2 <= 0:
            pares.append(None)
            continue
        pares.append(float(np.log(e1 / e2) / np.log(h1 / h2)))
    return pares


def ajuste_ordem(hs, erros):
    """Ajuste log E = log C + p log h nos pontos ainda decrescentes."""
    h = []
    e = []
    for hi, ei in zip(hs, erros):
        if ei is None or ei <= 0 or not np.isfinite(ei):
            continue
        if e and ei > e[-1] * 1.2:
            # A partir daqui o arredondamento ou a instabilidade dominam.
            break
        h.append(hi)
        e.append(ei)
    if len(h) < 2:
        return None
    coef = np.polyfit(np.log(h), np.log(e), 1)
    return {"p": float(coef[0]), "logC": float(coef[1]), "n_pontos": len(h)}


# ---------------------------------------------------------------------------
# Experimentos
# ---------------------------------------------------------------------------

def experimento_verificacao(L, C, V0) -> dict:
    """Confere as formas fechadas antes de usá-las como referência."""
    t = np.linspace(0.0, 5e-3, 20_001)
    relatorio = {}
    for nome, R in (("sub", 20.0), ("critico", 200.0), ("sobre", 500.0)):
        relatorio[nome] = {
            "R": R,
            "residuo": residuo_degrau_exato(t, V0, R, L, C),
            **{k: parametros_derivados(R, L, C)[k] for k in ("regime", "alpha", "Q", "zeta", "wd", "beta")},
        }
    info = parametros_derivados(20.0, L, C)
    for nome, w in (("ressonancia", info["w0"]), ("abaixo", 0.7 * info["w0"]), ("acima", 1.4 * info["w0"])):
        t_s = np.linspace(0.0, 8e-3, 80_001)
        amp, fase = fasor_corrente(V0, w, 20.0, L, C)
        y_fim = solucao_senoidal(np.array([8e-3]), V0, w, 20.0, L, C)[0, 0]
        permanente = amp * np.sin(w * 8e-3 + fase)
        relatorio[nome] = {
            "w": float(w),
            "residuo": residuo_senoidal(t_s, V0, w, 20.0, L, C),
            "amplitude": amp,
            "fase": fase,
            "erro_permanente_em_8ms": float(abs(y_fim - permanente)),
        }
    return relatorio


def experimento_convergencia(R, L, C, V0, tf, hs) -> dict:
    fonte = fonte_degrau(V0)
    f = fazer_campo(R, L, C, fonte)
    y0 = np.zeros(2)
    saida = {}
    for metodo, passo in PASSOS.items():
        serie = []
        for h in hs:
            t, y, nfe, ok, h_eff = integrar(passo, f, y0, 0.0, tf, h)
            ref = solucao_degrau(t, V0, R, L, C)
            met = erros_contra_referencia(y, ref)
            serie.append(
                {
                    "h_pedido": h,
                    "h": h_eff,
                    "passos": int(round(tf / h_eff)),
                    "avaliacoes": nfe,
                    "estavel": ok,
                    **met,
                }
            )
        erros = [item["max_i"] if item["estavel"] else None for item in serie]
        hs_eff = [item["h"] for item in serie]
        saida[metodo] = {
            "serie": serie,
            "ordens_consecutivas": ordem_observada(hs_eff, erros),
            "ajuste": ajuste_ordem(hs_eff, erros),
        }
    return saida


def experimento_estabilidade(casos, L, C, V0, tf, hs) -> dict:
    y0 = np.zeros(2)
    relatorio = {}
    for nome, R in casos:
        fonte = fonte_degrau(V0)
        f = fazer_campo(R, L, C, fonte)
        info = parametros_derivados(R, L, C)
        limites = {m: limite_estabilidade(m, R, L, C) for m in PASSOS}
        series = {}
        for metodo, passo in PASSOS.items():
            pontos = []
            for h in hs:
                t, y, nfe, ok, h_eff = integrar(passo, f, y0, 0.0, tf, h, limite=1e3)
                if ok and np.isfinite(y).all():
                    pico = float(np.nanmax(np.abs(y[:, 0])))
                else:
                    pico = None
                ref = solucao_degrau(t, V0, R, L, C)
                met = erros_contra_referencia(y, ref) if ok else {"max_i": None}
                pontos.append(
                    {
                        "h": h_eff,
                        "estavel_numerico": ok,
                        "pico_i": pico,
                        "max_i": met["max_i"],
                    }
                )
            series[metodo] = pontos
        relatorio[nome] = {"R": R, "info": info, "h_estabilidade": limites, "series": series}
    return relatorio


def experimento_regimes(L, C, V0, tf, h) -> dict:
    y0 = np.zeros(2)
    saida = {}
    for nome, R in (("subamortecido", 20.0), ("critico", 200.0), ("sobreamortecido", 500.0)):
        f = fazer_campo(R, L, C, fonte_degrau(V0))
        registro = {"R": R, "parametros": parametros_derivados(R, L, C), "metodos": {}}
        for metodo, passo in PASSOS.items():
            t, y, nfe, ok, h_eff = integrar(passo, f, y0, 0.0, tf, h)
            ref = solucao_degrau(t, V0, R, L, C)
            registro["metodos"][metodo] = {
                "h": h_eff,
                "estavel": ok,
                "avaliacoes": nfe,
                **erros_contra_referencia(y, ref),
            }
        saida[nome] = registro
    return saida


def experimento_senoidal(R, L, C, Vm, tf, h) -> dict:
    info = parametros_derivados(R, L, C)
    y0 = np.zeros(2)
    saida = {}
    for nome, fator in (("ressonancia", 1.0), ("abaixo", 0.7), ("acima", 1.4)):
        w = fator * info["w0"]
        f = fazer_campo(R, L, C, fonte_seno(Vm, w))
        amp, fase = fasor_corrente(Vm, w, R, L, C)
        item = {"w": float(w), "fator": fator, "amplitude": amp, "fase": fase, "metodos": {}}
        for metodo, passo in PASSOS.items():
            t, y, nfe, ok, h_eff = integrar(passo, f, y0, 0.0, tf, h)
            ref = solucao_senoidal(t, Vm, w, R, L, C)
            # Amplitude medida no último período.
            periodo = 2.0 * np.pi / w
            mascara = t >= tf - periodo
            amp_num = float(np.max(np.abs(y[mascara, 0]))) if ok else None
            item["metodos"][metodo] = {
                "h": h_eff,
                "estavel": ok,
                "avaliacoes": nfe,
                "amplitude_numerica": amp_num,
                "erro_amplitude": None if amp_num is None else abs(amp_num - amp),
                **erros_contra_referencia(y, ref),
            }
        saida[nome] = item
    return saida


def experimento_precisao(R, L, C, V0, tf, hs) -> dict:
    y0 = np.zeros(2)
    fonte = fonte_degrau(V0)
    saida = {}
    for nome, dtype in (("float64", np.float64), ("float32", np.float32)):
        f = fazer_campo(R, L, C, fonte, dtype=dtype)
        serie = []
        for h in hs:
            t, y, nfe, ok, h_eff = integrar(passo_rk4, f, y0, 0.0, tf, h, dtype=dtype)
            ref = solucao_degrau(t, V0, R, L, C)
            y64 = y.astype(float)
            met = erros_contra_referencia(y64, ref)
            serie.append({"h": h_eff, "estavel": ok, "avaliacoes": nfe, **met})
        erros = [item["max_i"] if item["estavel"] else None for item in serie]
        saida[nome] = {
            "serie": serie,
            "ajuste": ajuste_ordem([item["h"] for item in serie], erros),
        }
    return saida


def experimento_adaptativo(R, L, C, V0, tf) -> dict:
    f = fazer_campo(R, L, C, fonte_degrau(V0))
    y0 = np.zeros(2)
    tolerancias = [1e-4, 1e-6, 1e-8]
    resultados = []
    for rtol in tolerancias:
        t, y, hs, nfe, aceitos, rejeitados = rk4_adaptativo(
            f, y0, 0.0, tf, h0=2e-5, rtol=rtol, atol=rtol * 1e-3
        )
        ref = solucao_degrau(t, V0, R, L, C)
        met = erros_contra_referencia(y, ref)
        resultados.append(
            {
                "rtol": rtol,
                "atol": rtol * 1e-3,
                "avaliacoes": nfe,
                "aceitos": aceitos,
                "rejeitados": rejeitados,
                "h_min": float(np.min(hs)) if len(hs) else None,
                "h_mediano": float(np.median(hs)) if len(hs) else None,
                "h_max": float(np.max(hs)) if len(hs) else None,
                **met,
            }
        )
    return {"degrau": resultados}


def experimento_pulso(R, L, C) -> dict:
    """Pulso gaussiano curto: o passo adaptativo precisa cair e depois subir."""
    Vm, t0, sigma = 5.0, 1.0e-3, 40e-6
    tf = 4.0e-3
    fonte = fonte_pulso(Vm, t0, sigma)
    f = fazer_campo(R, L, C, fonte)
    y0 = np.zeros(2)
    # Referência: RK4 fixo com passo bem menor que a largura do pulso.
    t_ref, y_ref, nfe_ref, ok_ref, h_ref = integrar(passo_rk4, f, y0, 0.0, tf, 2.5e-7)
    t_grosso, y_grosso, _, ok_g, _ = integrar(passo_rk4, f, y0, 0.0, tf, 5.0e-7)
    # Concordância entre duas malhas finas, interpolando a mais grossa na fina.
    y_g_interp = np.column_stack(
        [np.interp(t_ref, t_grosso, y_grosso[:, 0]), np.interp(t_ref, t_grosso, y_grosso[:, 1])]
    )
    acordo = erros_contra_referencia(y_g_interp, y_ref)

    t_ad, y_ad, hs, nfe_ad, aceitos, rejeitados = rk4_adaptativo(
        f, y0, 0.0, tf, h0=5e-5, rtol=1e-6, atol=1e-8
    )
    y_ad_interp = np.column_stack(
        [np.interp(t_ref, t_ad, y_ad[:, 0]), np.interp(t_ref, t_ad, y_ad[:, 1])]
    )
    met = erros_contra_referencia(y_ad_interp, y_ref)

    # RK4 fixo com número de passos parecido com o adaptativo, para comparar custo.
    h_equiv = tf / max(aceitos, 1)
    t_fix, y_fix, nfe_fix, ok_fix, h_fix = integrar(passo_rk4, f, y0, 0.0, tf, h_equiv)
    y_fix_interp = np.column_stack(
        [np.interp(t_ref, t_fix, y_fix[:, 0]), np.interp(t_ref, t_fix, y_fix[:, 1])]
    )
    met_fix = erros_contra_referencia(y_fix_interp, y_ref)

    return {
        "Vm": Vm,
        "t0": t0,
        "sigma": sigma,
        "tf": tf,
        "referencia_h": h_ref,
        "referencia_avaliacoes": nfe_ref,
        "acordo_malhas_finas_max_i": acordo["max_i"],
        "adaptativo": {
            "rtol": 1e-6,
            "avaliacoes": nfe_ad,
            "aceitos": aceitos,
            "rejeitados": rejeitados,
            "h_min": float(np.min(hs)),
            "h_max": float(np.max(hs)),
            "h_mediano": float(np.median(hs)),
            **met,
        },
        "fixo_passo_equivalente": {
            "h": h_fix,
            "avaliacoes": nfe_fix,
            "estavel": ok_fix,
            **met_fix,
        },
        "_serie": {
            "t_ref": t_ref,
            "i_ref": y_ref[:, 0],
            "vc_ref": y_ref[:, 1],
            "t_ad": t_ad,
            "i_ad": y_ad[:, 0],
            "h_ad": hs,
            "vs": np.array([fonte(tt) for tt in t_ref]),
        },
    }


def custo_para_tolerancia(convergencia, tf, alvos) -> dict:
    """Estima, pelo ajuste E = C h^p, o custo para um erro maximo alvo."""
    tabela = {}
    for metodo, bloco in convergencia.items():
        ajuste = bloco["ajuste"]
        if not ajuste:
            continue
        p = ajuste["p"]
        C = float(np.exp(ajuste["logC"]))
        linhas = []
        for alvo in alvos:
            # C h^p = alvo => h = (alvo/C)^(1/p)
            if C <= 0 or p <= 0:
                continue
            h = (alvo / C) ** (1.0 / p)
            n = tf / h
            linhas.append(
                {
                    "erro_alvo": alvo,
                    "h": float(h),
                    "passos": float(n),
                    "avaliacoes": float(ESTAGIOS[metodo] * n),
                }
            )
        tabela[metodo] = {"C": C, "p": p, "linhas": linhas}
    return tabela


# ---------------------------------------------------------------------------
# Figuras
# ---------------------------------------------------------------------------

def figura_degrau(R, L, C, V0, tf, h, caminho: Path) -> None:
    f = fazer_campo(R, L, C, fonte_degrau(V0))
    t_fino = np.linspace(0.0, tf, 4000)
    ref = solucao_degrau(t_fino, V0, R, L, C)
    fig, eixos = plt.subplots(2, 1, figsize=(7.4, 6.2), sharex=True)
    eixos[0].plot(t_fino * 1e3, ref[:, 0] * 1e3, color=COR["analitica"], lw=2.0, label="Solução analítica")
    eixos[1].plot(t_fino * 1e3, ref[:, 1], color=COR["analitica"], lw=2.0, label="Solução analítica")
    for metodo, passo in PASSOS.items():
        t, y, _, ok, h_eff = integrar(passo, f, np.zeros(2), 0.0, tf, h)
        if not ok:
            continue
        eixos[0].plot(t * 1e3, y[:, 0] * 1e3, color=COR[metodo], lw=1.15, label=rotulo(metodo))
        eixos[1].plot(t * 1e3, y[:, 1], color=COR[metodo], lw=1.15, label=rotulo(metodo))
    eixos[0].set_ylabel("Corrente (mA)")
    eixos[1].set_ylabel("Tensão no capacitor (V)")
    eixos[1].set_xlabel("Tempo (ms)")
    eixos[0].set_title(f"Resposta ao degrau de 5 V, passo h = {h*1e6:.0f} µs")
    eixos[0].legend(loc="upper right", framealpha=0.92)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_erro_temporal(R, L, C, V0, tf, hs, caminho: Path) -> None:
    f = fazer_campo(R, L, C, fonte_degrau(V0))
    fig, eixos = plt.subplots(1, 2, figsize=(7.6, 3.6), sharey=False)
    for ax, h in zip(eixos, hs):
        for metodo, passo in PASSOS.items():
            t, y, _, ok, _ = integrar(passo, f, np.zeros(2), 0.0, tf, h)
            ref = solucao_degrau(t, V0, R, L, C)
            erro = np.abs(y[1:, 0] - ref[1:, 0])
            ax.semilogy(t[1:] * 1e3, np.maximum(erro, 1e-18), color=COR[metodo], lw=1.2, label=rotulo(metodo))
        ax.set_title(f"h = {h*1e6:.1f} µs")
        ax.set_xlabel("Tempo (ms)")
        ax.set_ylabel("Erro absoluto em i (A)")
        ax.legend(loc="lower right", framealpha=0.92)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_convergencia(convergencia, caminho: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    marcadores = {"euler": "o", "heun": "s", "rk4": "D"}
    for metodo, bloco in convergencia.items():
        hs, es = [], []
        for item in bloco["serie"]:
            if item["estavel"] and item["max_i"] is not None:
                hs.append(item["h"])
                es.append(item["max_i"])
        ax.loglog(hs, es, marker=marcadores[metodo], color=COR[metodo], lw=1.4, label=rotulo(metodo))
    # Retas de referência ancoradas no primeiro ponto do Euler, Heun e RK4.
    href = np.array([2.5e-7, 8e-6])
    for metodo, p, ls in (("euler", 1, ":"), ("heun", 2, "--"), ("rk4", 4, "-.")):
        serie = [it for it in convergencia[metodo]["serie"] if it["max_i"]]
        h0, e0 = serie[0]["h"], serie[0]["max_i"]
        ax.loglog(href, e0 * (href / h0) ** p, ls=ls, color=COR[metodo], alpha=0.55, lw=1.0)
    ax.set_xlabel("Passo h (s)")
    ax.set_ylabel("Erro máximo em i (A)")
    ax.set_title("Convergência no degrau subamortecido")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_estabilidade(estabilidade, caminho: Path) -> None:
    nomes = list(estabilidade.keys())
    fig, eixos = plt.subplots(1, len(nomes), figsize=(7.6, 3.8), sharey=True)
    if len(nomes) == 1:
        eixos = [eixos]
    for ax, nome in zip(eixos, nomes):
        bloco = estabilidade[nome]
        for metodo in PASSOS:
            hs, picos = [], []
            for item in bloco["series"][metodo]:
                hs.append(item["h"] * 1e6)
                pico = item["pico_i"]
                picos.append(np.nan if pico is None else min(pico, 50.0))
            ax.semilogy(hs, picos, marker="o", ms=3.5, color=COR[metodo], lw=1.2, label=rotulo(metodo))
            hcrit = bloco["h_estabilidade"][metodo] * 1e6
            ax.axvline(hcrit, color=COR[metodo], ls=":", lw=1.0, alpha=0.8)
        ax.set_title(f"R = {bloco['R']:.0f} \u03a9, Q = {bloco['info']['Q']:.1f}")
        ax.set_xlabel("Passo h (µs)")
        ax.set_ylabel("Pico de |i| no intervalo (A)")
        ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_regimes(L, C, V0, tf, h, caminho: Path) -> None:
    fig, eixos = plt.subplots(3, 1, figsize=(7.4, 7.4), sharex=True)
    casos = [
        ("Subamortecido, R = 20 \u03a9", 20.0),
        ("Crítico, R = 200 \u03a9", 200.0),
        ("Sobreamortecido, R = 500 \u03a9", 500.0),
    ]
    for ax, (titulo, R) in zip(eixos, casos):
        t = np.linspace(0.0, tf, 3000)
        ref = solucao_degrau(t, V0, R, L, C)
        f = fazer_campo(R, L, C, fonte_degrau(V0))
        ax.plot(t * 1e3, ref[:, 0] * 1e3, color=COR["analitica"], lw=2.0, label="Analítica")
        t_n, y_n, _, ok, _ = integrar(passo_rk4, f, np.zeros(2), 0.0, tf, h)
        if ok:
            ax.plot(t_n * 1e3, y_n[:, 0] * 1e3, color=COR["rk4"], lw=1.2, label="RK4")
        ax.set_ylabel("i (mA)")
        ax.set_title(titulo)
        ax.legend(loc="upper right")
    eixos[-1].set_xlabel("Tempo (ms)")
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_senoidal(R, L, C, Vm, tf, h, caminho: Path) -> None:
    info = parametros_derivados(R, L, C)
    fig, eixos = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
    casos = [
        ("Ressonância, \u03c9 = \u03c90", info["w0"]),
        ("Fora da ressonância, \u03c9 = 0,7 \u03c90", 0.7 * info["w0"]),
    ]
    for ax, (titulo, w) in zip(eixos, casos):
        t = np.linspace(0.0, tf, 5000)
        ref = solucao_senoidal(t, Vm, w, R, L, C)
        f = fazer_campo(R, L, C, fonte_seno(Vm, w))
        t_n, y, _, ok, _ = integrar(passo_rk4, f, np.zeros(2), 0.0, tf, h)
        ax.plot(t * 1e3, ref[:, 0] * 1e3, color=COR["analitica"], lw=1.8, label="Analítica")
        if ok:
            ax.plot(t_n * 1e3, y[:, 0] * 1e3, color=COR["rk4"], lw=1.05, label="RK4")
        amp, _ = fasor_corrente(Vm, w, R, L, C)
        ax.axhline(amp * 1e3, color="#7f8c8d", ls="--", lw=0.9, label="Amplitude permanente")
        ax.axhline(-amp * 1e3, color="#7f8c8d", ls="--", lw=0.9)
        ax.set_ylabel("Corrente (mA)")
        ax.set_title(titulo)
        ax.legend(loc="upper right")
    eixos[-1].set_xlabel("Tempo (ms)")
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_precisao(precisao, caminho: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for nome, cor, marca in (("float64", COR["rk4"], "D"), ("float32", COR["float32"], "o")):
        hs, es = [], []
        for item in precisao[nome]["serie"]:
            if item["max_i"] is not None and item["estavel"]:
                hs.append(item["h"])
                es.append(max(item["max_i"], 1e-18))
        ax.loglog(hs, es, marker=marca, color=cor, lw=1.4, label=nome)
    ax.set_xlabel("Passo h (s)")
    ax.set_ylabel("Erro máximo em i (A)")
    ax.set_title("RK4 em precisão dupla e simples")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def figura_pulso(pulso, caminho: Path) -> None:
    s = pulso["_serie"]
    fig, eixos = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=False)
    eixos[0].plot(s["t_ref"] * 1e3, s["vs"], color=COR["pulso"], lw=1.3, label="Pulso de excitação")
    eixos[0].plot(s["t_ref"] * 1e3, s["i_ref"] * 200, color=COR["analitica"], lw=1.3, label="Corrente de referência × 200")
    eixos[0].set_ylabel("Tensão (V) e corrente escalada")
    eixos[0].set_xlabel("Tempo (ms)")
    eixos[0].set_title("Resposta do tanque a um pulso curto")
    eixos[0].legend(loc="upper right")
    eixos[1].plot(s["t_ad"][:-1] * 1e3, s["h_ad"] * 1e6, color=COR["rk4"], lw=1.2)
    eixos[1].set_ylabel("Passo aceito h (µs)")
    eixos[1].set_xlabel("Tempo (ms)")
    eixos[1].set_title("Passo escolhido pelo RK4 adaptativo")
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)


def rotulo(metodo: str) -> str:
    return {"euler": "Euler", "heun": "Euler melhorado", "rk4": "RK4"}[metodo]


def limpar_serie(obj):
    """Remove vetores grandes antes de gravar o JSON."""
    if isinstance(obj, dict):
        return {k: limpar_serie(v) for k, v in obj.items() if not k.startswith("_")}
    if isinstance(obj, list):
        return [limpar_serie(v) for v in obj]
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def main() -> None:
    configurar_estilo()
    FIG.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)

    L = 10e-3
    C = 1e-6
    V0 = 5.0
    R = 20.0
    tf_conv = 5e-3
    tf_estab = 6e-3
    tf_seno = 8e-3
    hs_conv = [8e-6, 4e-6, 2e-6, 1e-6, 5e-7, 2.5e-7]
    hs_estab = [5e-6, 10e-6, 15e-6, 20e-6, 25e-6, 40e-6, 80e-6, 120e-6, 200e-6, 280e-6]
    hs_prec = [8e-6, 4e-6, 2e-6, 1e-6, 5e-7, 2.5e-7, 1.25e-7, 6.25e-8, 3.125e-8]

    print("Verificando soluções analíticas...")
    verificacao = experimento_verificacao(L, C, V0)
    for nome, item in verificacao.items():
        residuo = item["residuo"]
        print(f"  {nome:12s} resíduo = {residuo:.3e}")
        if residuo > 5e-3:
            raise SystemExit(f"Solução analítica inconsistente em {nome}: resíduo {residuo}")

    print("Convergência...")
    convergencia = experimento_convergencia(R, L, C, V0, tf_conv, hs_conv)
    print("Estabilidade...")
    estabilidade = experimento_estabilidade(
        [("Q5", 20.0), ("Q50", 2.0)], L, C, V0, tf_estab, hs_estab
    )
    print("Regimes de amortecimento...")
    regimes = experimento_regimes(L, C, V0, tf_conv, h=2e-6)
    print("Excitação senoidal...")
    senoidal = experimento_senoidal(R, L, C, V0, tf_seno, h=1e-6)
    print("Precisão float32/float64...")
    precisao = experimento_precisao(R, L, C, V0, tf=4e-3, hs=hs_prec)
    print("Passo adaptativo...")
    adaptativo = experimento_adaptativo(R, L, C, V0, tf_conv)
    print("Pulso...")
    pulso = experimento_pulso(R, L, C)
    custos = custo_para_tolerancia(convergencia, tf_conv, [1e-3, 1e-4, 1e-6, 1e-8])

    base = {
        "L_H": L,
        "C_F": C,
        "V0_V": V0,
        "R_principal_ohm": R,
        "tf_convergencia_s": tf_conv,
        "tf_estabilidade_s": tf_estab,
        "tf_senoidal_s": tf_seno,
        "circuito_principal": parametros_derivados(R, L, C),
        "circuito_Q50": parametros_derivados(2.0, L, C),
        "circuito_critico": parametros_derivados(200.0, L, C),
        "circuito_sobre": parametros_derivados(500.0, L, C),
        "epslon_float64": float(np.finfo(float).eps),
        "epslon_float32": float(np.finfo(np.float32).eps),
    }

    print("Gerando figuras...")
    figura_degrau(R, L, C, V0, tf_estab, 15e-6, FIG / "resposta_degrau.pdf")
    figura_erro_temporal(R, L, C, V0, tf_conv, [8e-6, 2e-6], FIG / "erro_temporal.pdf")
    figura_convergencia(convergencia, FIG / "convergencia.pdf")
    figura_estabilidade(estabilidade, FIG / "estabilidade.pdf")
    figura_regimes(L, C, V0, tf_conv, 2e-6, FIG / "regimes.pdf")
    figura_senoidal(R, L, C, V0, tf_seno, 1e-6, FIG / "senoidal.pdf")
    figura_precisao(precisao, FIG / "precisao.pdf")
    figura_pulso(pulso, FIG / "pulso.pdf")

    resultado = {
        "base": base,
        "verificacao": verificacao,
        "convergencia": convergencia,
        "estabilidade": estabilidade,
        "regimes": regimes,
        "senoidal": senoidal,
        "precisao": precisao,
        "adaptativo": adaptativo,
        "pulso": limpar_serie(pulso),
        "custos": custos,
    }
    destino = RES / "metricas.json"
    destino.write_text(json.dumps(limpar_serie(resultado), indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Métricas em {destino}")

    # Resumo legível no terminal.
    print("\n=== Limites de estabilidade (s) ===")
    for nome, bloco in estabilidade.items():
        print(nome, {k: f"{v:.6e}" for k, v in bloco["h_estabilidade"].items()})
    print("\n=== Ajuste de ordem ===")
    for metodo, bloco in convergencia.items():
        print(metodo, bloco["ajuste"], bloco["ordens_consecutivas"])
    print("\n=== Custos estimados ===")
    print(json.dumps(custos, indent=2))
    print("\n=== Adaptativo ===")
    print(json.dumps(adaptativo, indent=2))
    print("\n=== Pulso ===")
    print(json.dumps(limpar_serie(pulso), indent=2))
    print("\n=== Senoidal max_i ===")
    for nome, item in senoidal.items():
        print(nome, "amp", item["amplitude"])
        for m, met in item["metodos"].items():
            print("  ", m, met["max_i"], met["erro_amplitude"])
    print("\n=== Precisão (max_i) ===")
    for nome, bloco in precisao.items():
        print(nome, [(f"{s['h']:.3e}", s["max_i"]) for s in bloco["serie"]])


if __name__ == "__main__":
    main()
