"""Candidate fits of the logistic mixed model null, and the panels to try them on.

A is pyNei's: one explicit inverse of an individuals x individuals matrix per
linearization. B never forms an inverse while it iterates: it factors that
matrix with a Cholesky and solves, and the one quantity that needs more, the
trace of P times the kinship, comes from solving the kinship against that same
factorization, once for each step on tau instead of once per linearization.
Both stop at the same place, so their numbers should agree.
"""
import math, numpy
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.linalg.lapack import dtrtri
import pynei.gwas as g

TOL = g.GLMM_TOL
MAX_ITER = g.GLMM_MAX_ITER


def simulate(n_ind, n_var, seed=42, fst=0.1, family_size=4, num_pops=3, missing=0.0):
    """Families of full sibs inside subpops, as pyNei's reference panel is made."""
    rng = numpy.random.default_rng(seed)
    num_families = n_ind // family_size
    family_pops = rng.integers(0, num_pops, size=num_families)
    p_anc = rng.uniform(0.1, 0.9, size=n_var)
    a = p_anc * (1 - fst) / fst
    b = (1 - p_anc) * (1 - fst) / fst
    p_pop = numpy.stack([rng.beta(a, b) for _ in range(num_pops)])
    alleles = numpy.empty((n_var, n_ind, 2), dtype=numpy.int8)
    pops = numpy.repeat(family_pops, family_size)
    for fam, pop in enumerate(family_pops):
        parents = (rng.uniform(size=(2, n_var, 2)) < p_pop[pop][None, :, None]).astype(numpy.int8)
        for child in range(family_size):
            idx = fam * family_size + child
            for parent in range(2):
                picked = rng.integers(0, 2, size=n_var)
                alleles[:, idx, parent] = parents[parent, numpy.arange(n_var), picked]
    dosages = alleles.sum(axis=2).T.astype(float)
    p = dosages.mean(axis=0) / 2
    poly = (p > 0.05) & (p < 0.95)
    z = (dosages[:, poly] - 2 * p[poly]) / numpy.sqrt(2 * p[poly] * (1 - p[poly]))
    kinship = z @ z.T / z.shape[1]
    effects = rng.standard_normal(z.shape[1]) * math.sqrt(0.5 / z.shape[1])
    genetic = z @ effects
    cov1 = rng.standard_normal(n_ind)
    cov2 = rng.integers(0, 2, size=n_ind).astype(float)
    liability = (0.5 * cov1 + 0.8 * cov2 + 0.7 * pops + genetic
                 + rng.standard_normal(n_ind) * math.sqrt(0.5))
    y = (liability > numpy.quantile(liability, 0.6)).astype(float)
    design = numpy.column_stack([numpy.ones(n_ind), cov1, cov2])
    return y, design, kinship, dosages.T          # dosages: vars x individuals


# ---------------------------------------------------------------- A, pyNei's

def fit_a(y, design, kinship):
    null = g._GLMMNull(y, design, kinship)
    return {"tau": null.genetic_variance, "coefs": null.coefs,
            "projection": null.projection, "py": null.py}


# ------------------------------------------------- B, a Cholesky and no inverse

def _pql_for_tau_cholesky(y, design, kinship, tau, eta, mu, tol, max_iter, work):
    """The fixed effects and the linearization for one tau, without an inverse.

    It gives what the step on tau needs: the Cholesky of the covariance, the
    covariance solved against the design, and P applied to the working trait.
    """
    n = y.shape[0]
    for _ in range(max_iter):
        weights = mu * (1 - mu)
        working = eta + (y - mu) / weights
        sigma = tau * kinship
        sigma[numpy.diag_indices(n)] += 1 / weights
        chol = cho_factor(sigma, lower=True, overwrite_a=True, check_finite=False)
        work["choleskys"] += 1
        sid = cho_solve(chol, design, check_finite=False)          # n x c
        dsd = design.T @ sid
        coefs = numpy.linalg.solve(dsd, sid.T @ working)
        sw = cho_solve(chol, working, check_finite=False)
        pw = sw - sid @ numpy.linalg.solve(dsd, sid.T @ working)
        new_eta = design @ coefs + tau * (kinship @ pw)
        change = numpy.abs(new_eta - eta).max() / (numpy.abs(eta).max() + 1)
        eta = new_eta
        mu = g._expit(eta)
        if change < tol:
            break
    else:
        raise RuntimeError("The logistic mixed model did not converge")
    return coefs, eta, mu, chol, sid, dsd, pw, weights


def _apply_p(chol, sid, dsd, v):
    return cho_solve(chol, v, check_finite=False) - sid @ numpy.linalg.solve(dsd, sid.T @ v)


def _trace_p_kinship(chol, sid, dsd, kinship, work):
    """trace(P K) from the factorization: solve the kinship against it."""
    work["kinship_solves"] += 1
    sk = cho_solve(chol, kinship, check_finite=False)               # n x n
    first = numpy.trace(sk)
    second = numpy.trace(numpy.linalg.solve(dsd, sid.T @ (kinship @ sid)))
    return first - second


def fit_b(y, design, kinship, work=None):
    work = work if work is not None else {"choleskys": 0, "kinship_solves": 0}
    n, num_coefs = design.shape
    coefs, mu = g._fit_logistic(y, design)
    eta = design @ coefs
    weights = mu * (1 - mu)
    working = eta + (y - mu) / weights
    tau = float(numpy.var(working)) / 2
    too_small = too_big = None
    for _ in range(MAX_ITER):
        coefs, eta, mu, chol, sid, dsd, pw, weights = _pql_for_tau_cholesky(
            y, design, kinship, tau, eta, mu, TOL, MAX_ITER, work)
        kpw = kinship @ pw
        score = 0.5 * (pw @ kpw - _trace_p_kinship(chol, sid, dsd, kinship, work))
        ai = 0.5 * (kpw @ _apply_p(chol, sid, dsd, kpw))
        step = score / ai
        if abs(step) < TOL * (tau + TOL):
            break
        if score > 0:
            too_small = tau
        else:
            too_big = tau
        new_tau = tau + step
        if too_small is not None and too_big is not None:
            if not too_small < new_tau < too_big:
                new_tau = math.sqrt(too_small * too_big)
        elif new_tau <= 0:
            new_tau = tau / 4
        if new_tau < TOL:
            new_tau = 0.0
            if tau == 0.0:
                break
        tau = new_tau
    else:
        raise RuntimeError("The logistic mixed model did not converge")
    # the score test needs P itself, so it is formed once, at the end, from the
    # factorization of the last linearization, which is the one pyNei's P is of
    identity = numpy.eye(n)
    sigma_inv = cho_solve(chol, identity, check_finite=False)
    projection = g._projection(sigma_inv, design)
    return {"tau": tau, "coefs": coefs, "projection": projection, "py": y - mu,
            "work": work}


# ------------- C, B with the trace taken from the identity sigma = tau K + W^-1

def _trace_p_kinship_by_identity(chol, sid, dsd, kinship, tau, weights, work):
    """trace(P K), using that tau K = sigma - W^-1 makes the first term cheap.

    sigma^-1 K = (I - sigma^-1 W^-1) / tau, so the trace of sigma^-1 K is
    (n - trace(sigma^-1 W^-1)) / tau, and trace(sigma^-1 W^-1) is the squared
    length of the columns of L^-1 W^-1/2, one triangular solve where solving
    the kinship against the factorization takes two.
    """
    work["triangular_solves"] += 1
    n = kinship.shape[0]
    inv_l, info = dtrtri(chol[0], lower=1)
    if info:
        raise RuntimeError("the triangular factor could not be inverted")
    inv_l = numpy.tril(inv_l)
    first = (n - float(((inv_l * inv_l).sum(axis=0) / weights).sum())) / tau
    second = numpy.trace(numpy.linalg.solve(dsd, sid.T @ (kinship @ sid)))
    return first - second


def fit_c(y, design, kinship, work=None):
    work = work if work is not None else {"choleskys": 0, "triangular_solves": 0}
    n, num_coefs = design.shape
    coefs, mu = g._fit_logistic(y, design)
    eta = design @ coefs
    weights = mu * (1 - mu)
    working = eta + (y - mu) / weights
    tau = float(numpy.var(working)) / 2
    too_small = too_big = None
    for _ in range(MAX_ITER):
        coefs, eta, mu, chol, sid, dsd, pw, weights_used = _pql_for_tau_cholesky(
            y, design, kinship, tau, eta, mu, TOL, MAX_ITER, work)
        kpw = kinship @ pw
        trace = _trace_p_kinship_by_identity(chol, sid, dsd, kinship, tau, weights_used, work)
        score = 0.5 * (pw @ kpw - trace)
        ai = 0.5 * (kpw @ _apply_p(chol, sid, dsd, kpw))
        step = score / ai
        if abs(step) < TOL * (tau + TOL):
            break
        if score > 0:
            too_small = tau
        else:
            too_big = tau
        new_tau = tau + step
        if too_small is not None and too_big is not None:
            if not too_small < new_tau < too_big:
                new_tau = math.sqrt(too_small * too_big)
        elif new_tau <= 0:
            new_tau = tau / 4
        if new_tau < TOL:
            new_tau = 0.0
            if tau == 0.0:
                break
        tau = new_tau
    else:
        raise RuntimeError("The logistic mixed model did not converge")
    identity = numpy.eye(n)
    sigma_inv = cho_solve(chol, identity, check_finite=False)
    projection = g._projection(sigma_inv, design)
    return {"tau": tau, "coefs": coefs, "projection": projection, "py": y - mu,
            "work": work}


# ---------- D, C with the linearization run loosely while tau is still moving

def fit_d(y, design, kinship, work=None, loose=1e-3):
    """C, but the linearization is only converged tightly once tau has settled.

    The score that the step on tau is taken from is computed at the
    linearization, so a loose one early moves tau along a different path; the
    last steps use the tight tolerance, so the fit stops at the same place.
    """
    work = work if work is not None else {"choleskys": 0, "triangular_solves": 0}
    n, num_coefs = design.shape
    coefs, mu = g._fit_logistic(y, design)
    eta = design @ coefs
    weights = mu * (1 - mu)
    working = eta + (y - mu) / weights
    tau = float(numpy.var(working)) / 2
    too_small = too_big = None
    inner_tol = loose
    for _ in range(MAX_ITER):
        coefs, eta, mu, chol, sid, dsd, pw, weights_used = _pql_for_tau_cholesky(
            y, design, kinship, tau, eta, mu, inner_tol, MAX_ITER, work)
        kpw = kinship @ pw
        trace = _trace_p_kinship_by_identity(chol, sid, dsd, kinship, tau, weights_used, work)
        score = 0.5 * (pw @ kpw - trace)
        ai = 0.5 * (kpw @ _apply_p(chol, sid, dsd, kpw))
        step = score / ai
        if abs(step) < TOL * (tau + TOL) and inner_tol == TOL:
            break
        # once the step is small the linearization is converged tightly, and
        # from then on the fit is C's
        if abs(step) < loose * (tau + loose):
            inner_tol = TOL
        if score > 0:
            too_small = tau
        else:
            too_big = tau
        new_tau = tau + step
        if too_small is not None and too_big is not None:
            if not too_small < new_tau < too_big:
                new_tau = math.sqrt(too_small * too_big)
        elif new_tau <= 0:
            new_tau = tau / 4
        if new_tau < TOL:
            new_tau = 0.0
            if tau == 0.0:
                break
        tau = new_tau
    else:
        raise RuntimeError("The logistic mixed model did not converge")
    sigma_inv = cho_solve(chol, numpy.eye(n), check_finite=False)
    projection = g._projection(sigma_inv, design)
    return {"tau": tau, "coefs": coefs, "projection": projection, "py": y - mu,
            "work": work}
