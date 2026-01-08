import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def d32(bsd):
    """Compute d32 for a given 1-D array `bsd`.

    d32 = sum(x^3) / sum(x^2)
    """
    return np.sum(bsd**3) / np.sum(bsd**2)


def lognormal_pdf(x, mu, sigma):
    """Analytical lognormal PDF for x>0."""
    coef = 1.0 / (x * sigma * np.sqrt(2 * np.pi))
    exponent = -((np.log(x) - mu) ** 2) / (2 * sigma ** 2)
    return coef * np.exp(exponent)


def true_d32_from_lognormal(mu, sigma):
    """Closed-form population d32 for lognormal: E[X^3]/E[X^2]."""
    # E[X^k] = exp(k*mu + 0.5*k^2*sigma^2)
    return np.exp(mu + 2.5 * sigma ** 2)


if __name__ == "__main__":
    # Reproducible results
    rng = np.random.default_rng(42)
    # Printing / figure settings suitable for two-column A4 print
    column_width_mm = 86         # typical two-column width (mm)
    column_height_mm = 65        # chosen figure height (mm)
    column_width_in = column_width_mm / 25.4
    column_height_in = column_height_mm / 25.4
    fig_size = (column_width_in, column_height_in)
    dpi = 600                    # print-quality dpi
    save_figures = True         # set True to save high-quality PDFs
    
    png_path = Path("results/figures/png")
    svg_path = Path("results/figures/svg")
    eps_path = Path("results/figures/eps")

    # Configure matplotlib for print-quality figures
    fontfam = 'sans-serif'
    plt.rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "font.size": 10,
        # "axes.titlesize": 9,
        # "axes.labelsize": 8,
        # "legend.fontsize": 7,
        # "xtick.labelsize": 7,
        # "ytick.labelsize": 7,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        # "pdf.fonttype": 42,
        # "ps.fonttype": 42,
        "svg.fonttype": 'none',
        'font.family': fontfam,
        f'font.{fontfam}': ['arial', *plt.rcParams[f'font.{fontfam}']],
    })
    # Target a sharp peak (mode) around 2.0 with a small sigma
    mode_target = 2.0
    sigma = .7  # smaller sigma -> sharper peak
    mu = np.log(mode_target) + sigma ** 2  # solve mode = exp(mu - sigma^2)
    print(f"mu: {mu:.3f}")
    # Plot the theoretical PDF
    x = np.linspace(0.01, 20.0, 2000)
    pdf = lognormal_pdf(x, mu, sigma)

    # Choose increasing sample sizes to illustrate convergence
    sample_sizes = [20, 100, 500, 2000, 10000]
    max_n = max(sample_sizes)

    # Draw a single large pool and reuse prefixes for consistency
    samples_pool = rng.lognormal(mean=mu, sigma=sigma, size=max_n)

    # Figure 1: PDF + histograms for different sample sizes (print-sized)
    fig, ax = plt.subplots(1, 1, figsize=fig_size, layout="constrained")

    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_sizes)))
    bins = np.linspace(0, 20, 80)

    for n, c in zip(sample_sizes, colors):
        s = samples_pool[:n]
        ax.hist(s, bins=bins, density=True, color=c,
                label=f"n={n}", histtype='step')
    ax.plot(x, pdf, color="k", lw=1.5, label="Theoretical PDF")

    # plt.xlim(0, 6)
    ax.set_xlabel("Bubble diameter (mm)")
    ax.set_ylabel("Density")
    ax.set_title("Theoretical PDF and sample histograms")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.25)
    # ax.tight_layout()
    if save_figures:
        fig.savefig(eps_path / "figureD8_lognormal_hist.eps")
        fig.savefig(svg_path / "figureD8_lognormal_hist.svg")
        fig.savefig(png_path / "figureD8_lognormal_hist.png", dpi=dpi)
    # plt.show()

    # Figure 2: Evolution of d32 with increasing number of samples
    # Run multiple independent experiments (different seeds) to check robustness
    num_runs = 10
    seeds = rng.integers(0, 1_000_000_000, size=num_runs)
    print("Seeds used for runs:", seeds)

    n_range = np.arange(1, max_n + 1)
    d32_all = np.zeros((num_runs, max_n))
    median_all = np.zeros((num_runs, max_n))

    for i, seed in enumerate(seeds):
        r = np.random.default_rng(seed)
        samples = r.lognormal(mean=mu, sigma=sigma, size=max_n)
        print(f"Run {i+1}: max={samples.max():.1f}, mean={samples.mean():.2f}")
        cs2 = np.cumsum(samples ** 2)
        cs3 = np.cumsum(samples ** 3)
        d32_all[i, :] = cs3 / cs2
        # cumulative medians (prefix medians)
        # Note: this is O(n^2) in naive implementation but fine for these sizes
        median_all[i, :] = np.array([np.median(samples[:n]) for n in n_range])

    population_d32 = true_d32_from_lognormal(mu, sigma)
    population_median = np.exp(mu)

    # Combined figure: left = d32 evolution, right = median evolution
    double_fig_size = (fig_size[0], fig_size[1] * 2)
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=double_fig_size, constrained_layout=True)
    cmap = plt.cm.tab10

    # Left: d32 evolution
    for i in range(num_runs):
        ax0.plot(n_range, d32_all[i], lw=0.8, alpha=0.6, color=cmap(i % 10))

    d32_mean = d32_all.mean(axis=0)
    p10 = np.percentile(d32_all, 10, axis=0)
    p90 = np.percentile(d32_all, 90, axis=0)

    ax0.plot(n_range, d32_mean, color="k", lw=1.5, label="Mean $d_{32}$ across runs")
    ax0.fill_between(n_range, p10, p90, color="gray", alpha=0.25, label="10-90 percentile range")
    ax0.axhline(population_d32, color="k", lw=1.5, ls="--", label="Population $d_{32}$ = " + f"{population_d32:.4f}")

    # Highlight chosen sample sizes on the d32 plot
    # for n in sample_sizes:
    #     ax0.scatter(n, d32_mean[n - 1], s=20, color="C1")
    #     ax0.text(n, d32_mean[n - 1], f"  n={n}", va="center", fontsize=7)

    ax0.set_xscale("log")
    ax0.set_xlabel("Number of samples")
    ax0.set_ylabel("$d_{32}$ estimate")
    ax0.set_title("Evolution of $d_{32}$ estimate")
    ax0.grid(alpha=0.25)
    ax0.text(-.05, 1.05, "a", fontsize=9, fontweight='bold', transform=ax0.transAxes)
    # ax0.legend(frameon=False)

    # Right: median evolution
    for i in range(num_runs):
        ax1.plot(n_range, median_all[i], lw=0.8, alpha=0.6, color=cmap(i % 10))

    med_mean = median_all.mean(axis=0)
    med_p10 = np.percentile(median_all, 10, axis=0)
    med_p90 = np.percentile(median_all, 90, axis=0)

    ax1.plot(n_range, med_mean, color="k", lw=1.5, label="Mean median across runs")
    ax1.fill_between(n_range, med_p10, med_p90, color="gray", alpha=0.25, label="10-90 percentile range")
    ax1.axhline(population_median, color="k", lw=1.5, ls="--", label=f"Population median = {population_median:.4f}")

    # Highlight chosen sample sizes on the median plot
    # for n in sample_sizes:
    #     ax1.scatter(n, med_mean[n - 1], s=20, color="C1")
    #     ax1.text(n, med_mean[n - 1], f"  n={n}", va="center", fontsize=7)

    ax1.set_xscale("log")
    ax1.set_xlabel("Number of samples")
    ax1.set_ylabel("Median bubble size")
    ax1.set_title("Evolution of median bubble size")
    ax1.grid(alpha=0.25)
    ax1.text(-.05, 1.05, "b", fontsize=9, fontweight='bold', transform=ax1.transAxes)
    # ax1.legend(frameon=False)
    fig.legend(
        ax0.get_legend_handles_labels()[0],
        ["Mean $d_{32}$/$d_{median}$ across runs", "10-90 percentile range", "Population $d_{32}$/$d_{median}$"],
        # ax0.get_legend_handles_labels()[1],
        loc="outside lower center",
        # ncols=len(ax.get_legend_handles_labels()[0][1:])
    )

    if save_figures:
        fig.savefig(eps_path / "figureD9_d32_median_evolution.eps")
        fig.savefig(svg_path / "figureD9_d32_median_evolution.svg")
        fig.savefig(png_path / "figureD9_d32_median_evolution.png", dpi=dpi)
        
    # plt.show()
