// =====================================================================
// rlc_analise.sci
// Medidas de erro, ordem de convergência, limite de estabilidade e
// custo para um erro prescrito. Também grava os resultados em JSON.
// =====================================================================

function R = funcao_estab(metodo, z)
    // Função de estabilidade R(z) do método aplicado a y' = lambda y,
    // com z = h lambda. O passo é estável se |R(z)| <= 1.
    select metodo
    case "euler" then
        R = 1 + z;
    case "heun" then
        R = 1 + z + z.^2 / 2;
    case "rk4" then
        R = 1 + z + z.^2 / 2 + z.^3 / 6 + z.^4 / 24;
    end
endfunction

function hmax = limite_estab(metodo, p)
    // Maior passo para o qual os dois autovalores do circuito ficam dentro
    // da região de estabilidade absoluta. Varre h de 0,125 us a 5 ms e
    // para no primeiro h instável.
    hs = linspace(5e-3 / 40000, 5e-3, 40000);
    hmax = 0;
    for h = hs
        if or(abs(funcao_estab(metodo, h * p.lambda)) > 1 + 1e-9) then
            break
        end
        hmax = h;
    end
endfunction

function e = erro_max(y, ref)
    // Maior erro absoluto da corrente (e(1)) e da tensão (e(2)),
    // usando só os instantes em que as duas soluções são finitas.
    ok = and(~isnan(y), "r") & and(~isnan(ref), "r");
    if sum(ok) < 3 then
        e = [%nan, %nan];
        return
    end
    e = [max(abs(y(1, ok) - ref(1, ok))), max(abs(y(2, ok) - ref(2, ok)))];
endfunction

function q = ordem_observada(hs, es)
    // Ordem entre passos vizinhos: log(E1/E2) / log(h1/h2).
    q = [];
    for k = 1:size(hs, "*") - 1
        q($ + 1) = log(es(k) / es(k + 1)) / log(hs(k) / hs(k + 1));
    end
    q = q';
endfunction

function aj = ajuste_ordem(hs, es)
    // Reta de mínimos quadrados log E = log C + p log h.
    // Para quando o erro volta a crescer (arredondamento ou instabilidade).
    h = [];
    e = [];
    for k = 1:size(hs, "*")
        if isnan(es(k)) | es(k) <= 0 then
            continue
        end
        if size(e, "*") > 0 then
            if es(k) > e($) * 1.2 then
                break
            end
        end
        h($ + 1) = hs(k);
        e($ + 1) = es(k);
    end
    M = [log(h), ones(h)];
    c = M \ log(e);
    aj = struct("p", c(1), "logC", c(2), "n_pontos", size(h, "*"));
endfunction

function c = custo_para_erro(aj, estagios, tf, alvos)
    // Pelo ajuste E = C h^p, o passo para um erro alvo é h = (E/C)^(1/p).
    // O custo é o número de avaliações de f: estágios x (tf / h).
    C = exp(aj.logC);
    h = (alvos / C) .^ (1 / aj.p);
    passos = tf ./ h;
    c = struct("erro_alvo", alvos, "h", h, ..
               "passos", passos, "avaliacoes", estagios * passos);
endfunction

// ---------------------------------------------------------------------
// Gravação em JSON (feita à mão para funcionar em versões antigas do
// Scilab, que não têm toJSON).
// ---------------------------------------------------------------------

function s = json_numero(x)
    if isnan(x) | isinf(x) then
        s = "null";
    else
        s = msprintf("%.17g", x);
    end
endfunction

function s = json_valor(v, recuo)
    // Converte struct, list, número, vetor, texto ou booleano em JSON.
    nl = ascii(10);
    dentro = recuo + "  ";
    select typeof(v)
    case "constant" then
        if size(v, "*") == 1 then
            s = json_numero(v);
        else
            partes = [];
            for x = v(:)'
                partes($ + 1) = json_numero(x);
            end
            s = "[" + strcat(partes, ", ") + "]";
        end
    case "boolean" then
        if v then
            s = "true";
        else
            s = "false";
        end
    case "string" then
        s = """" + v + """";
    case "st" then
        campos = fieldnames(v);
        partes = [];
        for k = 1:size(campos, "*")
            partes($ + 1) = dentro + """" + campos(k) + """: " ..
                            + json_valor(v(campos(k)), dentro);
        end
        s = "{" + nl + strcat(partes, "," + nl) + nl + recuo + "}";
    case "list" then
        partes = [];
        for k = 1:size(v)
            partes($ + 1) = dentro + json_valor(v(k), dentro);
        end
        s = "[" + nl + strcat(partes, "," + nl) + nl + recuo + "]";
    else
        s = "null";
    end
endfunction

function gravar_json(v, arquivo)
    mputl(json_valor(v, ""), arquivo);
endfunction
