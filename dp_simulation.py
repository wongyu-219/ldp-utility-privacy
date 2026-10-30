"""
Section 5.4 experiment for:
  "Design of a User-Behavior Data Collection and Analysis System Applying
   Differential Privacy" (Won Gyu Lee, Hanyang University, June 2026)

Main experiment: Randomized Response (Scheme A: answer truthfully with
probability p, otherwise flip a fair coin), estimating a population
proportion pi = 0.30 under epsilon in {0.1, 0.5, 1.0, 2.0} and
n in {1,000, 10,000, 100,000}, 100 repetitions per condition.

Scheme A relations (see the erratum note for section 4.1):
  P(report=1 | true=1) = (1+p)/2,  P(report=1 | true=0) = (1-p)/2
  epsilon = ln((1+p)/(1-p))   <=>   p = (e^eps - 1)/(e^eps + 1)
  estimator: pi_hat = (y_hat - (1-p)/2) / p          (as in the report)

Supplementary experiment: Laplace mechanism for a mean on [0,1]
(section 4.2), true mean ~= 0.3, same epsilon-n grid, to check
SE ~ (b-a)/(eps*sqrt(n)) up to its constant.

Dependencies: numpy, pandas, matplotlib. Single fixed seed for
reproducibility.
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration (exactly section 5.4)
# ----------------------------------------------------------------------
SEED = 20260604
PI = 0.30                     # true proportion
N_POP = 1_000_000             # synthetic population size
K_POP = int(round(PI * N_POP))  # exactly 300,000 users with attribute = 1
EPSILONS = [0.1, 0.5, 1.0, 2.0]
NS = [1_000, 10_000, 100_000]
R = 100                       # repetitions per condition
Z = 1.959963984540054         # 97.5% normal quantile

rng = np.random.default_rng(SEED)

# Okabe-Ito colorblind-safe palette, fixed order by ascending epsilon
PALETTE = {0.1: "#0072B2", 0.5: "#E69F00", 1.0: "#009E73", 2.0: "#D55E00"}
MARKERS = {0.1: "o", 0.5: "s", 1.0: "^", 2.0: "D"}


def p_from_eps(eps: float) -> float:
    """Scheme A: p = (e^eps - 1)/(e^eps + 1) = tanh(eps/2)."""
    return (np.exp(eps) - 1.0) / (np.exp(eps) + 1.0)


# ----------------------------------------------------------------------
# Main experiment: Randomized Response (Scheme A)
# ----------------------------------------------------------------------
def run_rr() -> pd.DataFrame:
    rows = []
    for eps in EPSILONS:
        p = p_from_eps(eps)
        q1 = (1.0 + p) / 2.0          # P(report=1 | true=1)
        q0 = (1.0 - p) / 2.0          # P(report=1 | true=0)
        qbar = PI * q1 + (1.0 - PI) * q0
        for n in NS:
            est = np.empty(R)
            abs_err = np.empty(R)
            ci_width = np.empty(R)
            covered = np.zeros(R, dtype=bool)
            for r in range(R):
                # draw n users without replacement from the fixed population:
                # the number of true-1 users in the sample is hypergeometric
                k = rng.hypergeometric(K_POP, N_POP - K_POP, n)
                yes = rng.binomial(k, q1) + rng.binomial(n - k, q0)
                y_hat = yes / n
                pi_hat = (y_hat - (1.0 - p) / 2.0) / p   # report's estimator
                # 95% CI by normal approximation:
                # Var(y_hat) ~= y_hat(1-y_hat)/n, Var(pi_hat) = Var(y_hat)/p^2
                se_hat = np.sqrt(max(y_hat * (1.0 - y_hat), 0.0) / n) / p
                half = Z * se_hat
                est[r] = pi_hat
                abs_err[r] = abs(pi_hat - PI)
                ci_width[r] = 2.0 * half
                covered[r] = (pi_hat - half) <= PI <= (pi_hat + half)

            # theoretical SEs (two variants requested):
            # (a) conditional mixture only (ignores sampling of true values)
            se_mix = np.sqrt(
                (q1 * (1 - q1) * PI + q0 * (1 - q0) * (1 - PI)) / n
            ) / p
            # (b) total via qbar(1-qbar)/(n p^2) -- includes the binomial
            #     variance of the sample's true composition
            se_tot = np.sqrt(qbar * (1 - qbar) / n) / p

            rows.append(
                dict(
                    epsilon=eps,
                    p=p,
                    n=n,
                    mae=abs_err.mean(),
                    median_abs_err=np.median(abs_err),
                    p90_abs_err=np.quantile(abs_err, 0.90),
                    sd_empirical=est.std(ddof=1),
                    se_theory_mixture=se_mix,
                    se_theory_total=se_tot,
                    mean_ci_width=ci_width.mean(),
                    coverage=covered.mean(),
                    pct_outside_01=float(np.mean((est < 0) | (est > 1))),
                    mean_estimate=est.mean(),
                )
            )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Supplementary: Laplace mechanism for a bounded mean (section 4.2)
# ----------------------------------------------------------------------
def run_laplace() -> pd.DataFrame:
    # fixed synthetic population of per-user values on [0,1], mean ~ 0.3
    pop = rng.beta(3.0, 7.0, size=N_POP)   # mean 0.3 by construction
    true_mean = pop.mean()
    var_pop = pop.var()
    rows = []
    for eps in EPSILONS:
        b_scale = 1.0 / eps                 # Delta f = (b - a) = 1
        for n in NS:
            abs_err = np.empty(R)
            est = np.empty(R)
            for r in range(R):
                idx = rng.choice(N_POP, size=n, replace=False)
                noisy = pop[idx] + rng.laplace(0.0, b_scale, size=n)
                m = noisy.mean()
                est[r] = m
                abs_err[r] = abs(m - true_mean)
            rows.append(
                dict(
                    epsilon=eps,
                    n=n,
                    true_mean=true_mean,
                    mae=abs_err.mean(),
                    sd_empirical=est.std(ddof=1),
                    se_report=1.0 / (eps * np.sqrt(n)),          # (b-a)/(eps sqrt(n))
                    se_noise_exact=np.sqrt(2.0) / (eps * np.sqrt(n)),
                    se_total_theory=np.sqrt(var_pop / n + 2.0 / (eps**2 * n)),
                )
            )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def style_axes(ax):
    ax.grid(True, which="major", color="#d9d9d9", linewidth=0.6, alpha=0.8)
    ax.grid(True, which="minor", color="#ededed", linewidth=0.4, alpha=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#9a9a9a")
    ax.tick_params(colors="#444444", labelsize=9)


def fig1_mae_vs_n(df: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=200)
    for eps in EPSILONS:
        sub = df[df.epsilon == eps].sort_values("n")
        ax.plot(
            sub.n, sub.mae,
            color=PALETTE[eps], marker=MARKERS[eps], markersize=5,
            linewidth=1.8, label=f"ε = {eps}",
        )
        ax.annotate(
            f"ε={eps}", xy=(sub.n.iloc[-1], sub.mae.iloc[-1]),
            xytext=(6, 0), textcoords="offset points",
            fontsize=8.5, color="#333333", va="center",
        )
    # slope -1/2 reference line, anchored to eps=1.0 at n=1,000
    ref = df[(df.epsilon == 1.0) & (df.n == 1_000)].mae.iloc[0] * 1.8
    ns = np.array(NS, dtype=float)
    ax.plot(
        ns, ref * (ns / ns[0]) ** -0.5,
        color="#777777", linestyle="--", linewidth=1.2,
        label="slope −1/2 reference",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sample size n (log scale)", fontsize=10, color="#222222")
    ax.set_ylabel("MAE of π̂ (log scale)", fontsize=10, color="#222222")
    ax.set_title(
        "Figure 1. MAE vs. sample size (randomized response,\nπ = 0.30, 100 runs per condition)",
        fontsize=10, color="#222222", loc="left",
    )
    style_axes(ax)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.margins(x=0.08)
    fig.tight_layout()
    fig.subplots_adjust(right=0.88)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def fig2_mae_vs_eps(df: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=200)
    n_colors = {1_000: "#0072B2", 10_000: "#E69F00", 100_000: "#009E73"}
    n_markers = {1_000: "o", 10_000: "s", 100_000: "^"}
    for n in NS:
        sub = df[df.n == n].sort_values("epsilon")
        ax.plot(
            sub.epsilon, sub.mae,
            color=n_colors[n], marker=n_markers[n], markersize=5,
            linewidth=1.8, label=f"n = {n:,}",
        )
        ax.annotate(
            f"n={n:,}", xy=(sub.epsilon.iloc[-1], sub.mae.iloc[-1]),
            xytext=(6, 0), textcoords="offset points",
            fontsize=8.5, color="#333333", va="center",
        )
    # 1/eps reference, anchored to n=10,000 at eps=0.5
    ref = df[(df.epsilon == 0.5) & (df.n == 10_000)].mae.iloc[0] * 1.8
    eps_grid = np.array(EPSILONS)
    ax.plot(
        eps_grid, ref * (eps_grid / 0.5) ** -1.0,
        color="#777777", linestyle="--", linewidth=1.2, label="1/ε reference",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(EPSILONS)
    ax.set_xticklabels([str(e) for e in EPSILONS])
    ax.set_xlabel("Privacy budget ε (log scale)", fontsize=10, color="#222222")
    ax.set_ylabel("MAE of π̂ (log scale)", fontsize=10, color="#222222")
    ax.set_title(
        "Figure 2. MAE vs. privacy budget (randomized response,\nπ = 0.30, 100 runs per condition)",
        fontsize=10, color="#222222", loc="left",
    )
    style_axes(ax)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.margins(x=0.08)
    fig.tight_layout()
    fig.subplots_adjust(right=0.86)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
if __name__ == "__main__":
    import os
    os.makedirs("results", exist_ok=True)
    rr = run_rr()
    rr.to_csv("results/results.csv", index=False)
    lap = run_laplace()
    lap.to_csv("results/laplace_results.csv", index=False)
    fig1_mae_vs_n(rr, "results/fig1_mae_vs_n.png")
    fig2_mae_vs_eps(rr, "results/fig2_mae_vs_eps.png")

    pd.set_option("display.width", 200)
    pd.set_option("display.float_format", lambda v: f"{v:.5f}")
    print("=== Randomized Response (Scheme A) ===")
    print(rr.to_string(index=False))
    print("\n=== Laplace mechanism (supplementary) ===")
    print(lap.to_string(index=False))
