## functions for fitting models to estimated statistics

from datetime import datetime
import demes
import numpy as np
import moments
from moments.Demes import Inference
import scipy
import sys
import os
import pickle

from . import utils
from .datastrucs import DplusStats


# TODO write own getter/setter functions for fitting u


_out_of_bounds = 1e10
_counter = 0


def load_statistics(data_file, graph_file):

    data = pickle.load(open(data_file, "rb+"))
    pop_ids = data["pop_ids"]
    to_pops = graph_data_overlap(graph_file, pop_ids)
    bins = data["bins"]
    means = utils.subset_means(data["means"], pop_ids, to_pops)
    varcovs = utils.subset_varcovs(data["varcovs"], pop_ids, to_pops)

    return to_pops, bins, means, varcovs


def compute_bin_stats(
    graph,
    sampled_demes=None, 
    sample_times=None, 
    u=None,
    bins=None,
    approx="simpsons",
    phased=False
):
    """
    From a Demes graph, compute expected D+ in bins using moments.LD and a
    given approximation method. 

    :returns: A DPlusStats instance holding computed statistics.
    """
    if approx not in ("midpoint", "trapezoid", "simpsons"):
        raise ValueError("unrecognized approximation method")
    
    if isinstance(graph, str):
        graph = demes.load(graph)

    if approx == "midpoint":
        midpoints = (bins[:-1] + bins[1:]) / 2
        model = DplusStats.from_moments(
            graph, 
            sampled_demes, 
            sample_times=sample_times, 
            rs=midpoints,
            u=u,
            phased=phased
        )

    elif approx == "trapezoid":
        raise ValueError("this method is too inaccurate")
        y_edges = DplusStats.from_moments(
            graph, 
            sampled_demes, 
            sample_times=sample_times, 
            rs=bins, 
            u=u,
            phased=phased
        )
        y = [(y0 + y1) / 2 for y0, y1 in zip(y_edges[:-2], y_edges[1:-1])]
        y.append(y_edges[-1])
        model = DplusStats(y, pop_ids=sampled_demes)

    elif approx == "simpsons":
        y_edges = DplusStats.from_moments(
            graph, 
            sampled_demes, 
            sample_times=sample_times, 
            rs=bins, 
            u=u,
            phased=phased
        )
        midpoints = (bins[:-1] + bins[1:]) / 2
        y_mids = DplusStats.from_moments(
            graph, 
            sampled_demes, 
            sample_times=sample_times,
            rs=midpoints, 
            u=u,
            phased=phased
        )        
        y = [
            (y_edges[i] + 4 * y_mids[i] + y_edges[i + 1]) / 6 
            for i in range(len(midpoints))
        ]
        y.append(y_edges[-1])
        model = DplusStats(y, pop_ids=sampled_demes)

    else:
        raise ValueError("invalid approximation method!")

    return model


## optimization functions


def _object_func(
    params,
    builder,
    options,
    means,
    varcovs,
    sampled_demes=None,
    sample_times=None,
    u=None,
    bins=None,
    lower_bounds=None,
    upper_bounds=None,
    constraints=None,
    verbose=None,
    one_locus=False,
    use_afs=False,
    afs=None,
    stream=sys.stdout
):
    
    if lower_bounds is not None and np.any(params < lower_bounds):
        return -_out_of_bounds
    elif upper_bounds is not None and np.any(params > upper_bounds):
        return -_out_of_bounds
    elif constraints is not None and np.any(constraints(params) <= 0):
        return -_out_of_bounds

    global _counter
    _counter += 1    

    builder = Inference._update_builder(builder, options, params)
    graph = demes.Graph.fromdict(builder)
    model = compute_bin_stats(
        graph, 
        sampled_demes,
        sample_times=sample_times,
        u=u,
        bins=bins,
        phased=False
    )
    ll = composite_ll(model, means, varcovs, one_locus=one_locus)

    if use_afs:
        ll_afs = 0

        ll += ll_afs
        raise ValueError("not implemented")
    
    if verbose > 0 and _counter % verbose == 0:
        param_str = "".join(["{:>11}".format(p) for p in format_params(params)])
        report_str = "{n:<5g}{l:>10g}  [{p}]".format(
            n=_counter, l=np.round(ll, 2), p=param_str
        )
        print(report_str)

    return -ll


def _object_func_log(log_p, *args, **kwargs):
    
    p = np.exp(log_p - 1)
    return _object_func(p, *args, **kwargs)


def optimize(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=None,
    bins=None,
    u=None,
    method="fmin",
    max_iter=1000,
    max_calls=None,
    log=False,
    verbose=1,
    overwrite=False,
    output=None,
    one_locus=False,
    perturb=False,
    stream=sys.stdout
):
    """

    """
    builder = Inference._get_demes_dict(graph_file)
    options = Inference._get_params_dict(param_file)
    params_bounds = Inference._set_up_params_and_bounds(options, builder)
    param_names, params_0, lower_bounds, upper_bounds = params_bounds
    constraints = Inference._set_up_constraints(options, param_names)

    if u is None:
        raise ValueError("you must provide `u`")
    if pop_ids is None:
        raise ValueError("you must provide `pop_ids`")
    
    if perturb > 0: 
        params_0 = Inference._perturb_params_constrained(
            params_0, 
            perturb, 
            lower_bound=lower_bounds, 
            upper_bound=upper_bounds,
            cons=constraints
        )
    if log:
        objective = _object_func_log
        params_0 = np.log(params_0) + 1
    else:
        objective = _object_func
    
    print(get_time(), f"Fitting D+ to data for {pop_ids}")
    param_name_str = "".join(["{:>11}".format(p) for p in param_names])
    init_param_str = "{n:<5}{l:>10}  [{p}]".format(
        n="", l="Params", p=param_name_str
    )
    print(init_param_str)
    param_0_str = "".join(["{:>11}".format(p) for p in format_params(params_0)])
    init_param_str = "{n:<5}{l:>10}  [{p}]".format(
        n="Calls", l="LL", p=param_0_str
    )
    print(init_param_str)

    deme_names = [d["name"] for d in builder["demes"]]
    sampled_demes = [] 
    sample_times = []
    for pop in pop_ids: 
        assert pop in deme_names
        idx = deme_names.index(pop)
        sample_times.append(builder["demes"][idx]["epochs"][-1]["end_time"])
        sampled_demes.append(pop)
    
    warn = None
    args = (
        builder,
        options,
        means,
        varcovs,
        sampled_demes,
        sample_times,
        u,
        bins,
        lower_bounds,
        upper_bounds,
        constraints,
        verbose,
        one_locus,
        False,
        None,
        stream
    )
    
    methods = ['fmin', 'powell', 'bfgs', 'lbfgsb']
    if method not in methods:
        raise ValueError(f"{method} is not a valid method")
    
    if method == 'fmin':
        ret = scipy.optimize.fmin(
            objective,
            params_0,
            args=args,
            maxiter=max_iter,
            maxfun=max_calls,
            full_output=True
        )
        fit_params, fopt, num_iter, func_calls, flag = ret[:5]

    elif method == 'powell':
        ret = scipy.optimize.fmin_powell(
            objective,
            params_0,
            args=args,
            maxiter=max_iter,
            full_output=True,
        )
        fit_params, fopt, _, num_iter, func_calls, flag = ret[:6]

    elif method == 'bfgs':
        if log:
            epsilon = 1e-3
        else:
            epsilon = None
        ret = scipy.optimize.fmin_bfgs(
            objective,
            params_0,
            args=args,
            maxiter=max_iter,
            epsilon=epsilon,
            disp=False,
            full_output=True
        )
        fit_params, fopt, _, __, func_calls, grad_calls, flag = ret[:7]
        num_iter = grad_calls

    elif method == 'lbfgsb':
        if log:
            bounds = list(
                zip(np.log(lower_bounds) + 1, np.log(upper_bounds) + 1)
            )
            epsilon = 1e-5
            pgtol = 1e-5
        else:
            bounds = list(zip(lower_bounds, upper_bounds))
            epsilon = 1e-8
            pgtol = 1e-5
        ret = scipy.optimize.fmin_l_bfgs_b(
            objective,
            params_0,
            args=args,
            maxiter=max_iter,
            bounds=bounds,
            epsilon=epsilon,
            pgtol=pgtol,
            approx_grad=True
        )
        fit_params, fopt, output_dict = ret
        num_iter = output_dict['nit']
        flag = output_dict['warnflag']
        warn = output_dict["task"]

    else:
        return

    if log: 
        fit_params = np.exp(fit_params - 1)

    ll = -fopt

    print(f"Log-likelihood:\t{ll:.3}")
    print("Fitted parameters:")
    for name, value in zip(param_names, fit_params):
        print(f"{name}\t{value:.3}")

    global _counter

    if output is not None:
        builder = Inference._update_builder(builder, options, fit_params)
        graph = demes.Graph.fromdict(builder)
        if overwrite is False and os.path.isfile(output):
            print(f"{output} already exists: printing model")
            print(str(graph))
        else:
            info = dict(
                method=method,
                objective_func=objective.__name__,
                fopt=-ll,
                max_iter=max_iter,
                num_iter=num_iter,
                flag=flag,
                warn=warn,
                u=u
            )
            graph.metadata['opt_info'] = info
            demes.dump(graph, output)

    _counter = 0
    
    return param_names, fit_params, ll


def print_status(n_calls, ll, params):
    """
    Print the number of function calls, the log-likelihood, and the current 
    parameter values.
    """
    t = utils.get_time()
    _n = f'{n_calls:<4}'
    if isinstance(ll, float):
        _ll = f'{np.round(ll, 2):>10}'
    else:
        _ll = f'{ll:>10}'
    fmt_p = []
    for x in params:
        if isinstance(x, str):
            fmt_p.append(f'{x:>10}')
        else:
            if x > 1:
                fmt_p.append(f'{np.round(x, 1):>10}')
            else:
                sci = np.format_float_scientific(x, 2, trim='k')
                fmt_p.append(f'{sci:>10}')
    _p = ''.join(fmt_p)
    print(t, _n, _ll, '[', _p, ']')

    return


def format_params(params):
    """
    returns strings
    
    """
    formatted = []
    for param in params:
        if param >= 1:
            formatted.append(str(np.round(param, 1)))
        elif param >= 1e-3:
            formatted.append(np.format_float_positional(param, precision=3))
        else:
            formatted.append(np.format_float_scientific(param, precision=2))

    return formatted


## computing log-likelihoods


_inv_varcov_cache = dict()


def composite_ll(model, means, varcovs, one_locus=False):
    """
    Compute the sum of log-likelihoods across bins.
    """
    if one_locus:
        ll = ll_per_bin(model, means, varcovs).sum()
    else:
        ll = ll_per_bin(model[:-1], means[:-1], varcovs[:-1]).sum()

    return ll


def ll_per_bin(xs, mus, varcovs):
    """
    Compute LL in each bin and return array of bin LLs
    """
    n_bins = len(xs)
    if len(mus) != n_bins or len(varcovs) != n_bins:
        raise ValueError("Data, model and varcovs must have the same length")
    bin_ll = np.zeros(n_bins, dtype=np.float64)
    for i in range(n_bins):
        if (
            i in _inv_varcov_cache  
            and np.all(_inv_varcov_cache[i]["varcov"] == varcovs[i])
        ):
            inv_varcov = _inv_varcov_cache[i]["inv_varcov"]
        else:
            inv_varcov = np.linalg.inv(varcovs[i])
            add_to_cache = {"varcov": varcovs[i], "inv_varcov": inv_varcov}
            _inv_varcov_cache[i] = add_to_cache
        bin_ll[i] = _ll(xs[i], mus[i], inv_varcov)

    return bin_ll


def _ll(x, mu, inv_varcov):
    """
    Compute the log of the multivariate gaussian function with means `mu`, 
    pre-inverted covariance matrix `inv_cov` at `x`.

    :param x: Array at which to evaluate the function.
    :type x: np.ndarray, shape (n,)
    :param mu: Array specifying mean parameters for the distribution.
    :type mu: np.ndarray, shape (n,)
    :inv_varcov: Pre-inverted covariance matrix parameterizing the distribution
    :type inv_varcov: np.ndarray, shape (n, n)

    :returns: Log of the multivariate gaussian law.
    """
    return -1.0 / 2.0 * np.matmul(np.matmul(x - mu, inv_varcov), x - mu)


def _log_multivariate_normal_pdf(x, mu, varcov):

    f = _multivariate_normal_pdf(x, mu, varcov)
    return np.log(f)


def _multivariate_normal_pdf(x, mu, varcov):

    k = len(x)
    inv_varcov = np.linalg.inv(varcov)
    f = (
        np.exp(-1.0 / 2.0 * np.matmul(np.matmul(x - mu, inv_varcov), x - mu)) 
        / np.sqrt(np.linalg.det(varcov) * (2.0 * np.pi) ** k) 
    )
    return f


## Computing uncertainties; estimating standard errors


def compute_uncerts(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=None,
    bins=None,
    u=None,
    bootstrap_reps=None,
    delta=0.01,
    method="godambe"
):
    """
    
    """
    if method not in ("godambe", "fisher"):
        raise ValueError("invalid method")

    builder = Inference._get_demes_dict(graph_file)
    options = Inference._get_params_dict(param_file)
    params_bounds = Inference._set_up_params_and_bounds(options, builder)
    param_names, params_0, lower_bounds, upper_bounds = params_bounds
    deme_names = [d["name"] for d in builder["demes"]]
    sampled_demes = [] 
    sample_times = []
    for pop in pop_ids: 
        assert pop in deme_names
        idx = deme_names.index(pop)
        sample_times.append(builder["demes"][idx]["epochs"][-1]["end_time"])
        sampled_demes.append(pop)

    model_args = (
        builder, 
        options, 
        sampled_demes, 
        sample_times, 
        bins, 
        u
    )

    def model_func(params, args=()):
        # 
        builder, options, sampled_demes, sample_times, bins, u = args
        builder = Inference._update_builder(builder, options, params)
        graph = demes.Graph.fromdict(builder)
        model = compute_bin_stats(
            graph, 
            sampled_demes,
            sample_times=sample_times,
            u=u,
            bins=bins,
            phased=False
        )

        return model

    if method == "fisher":
        H = compute_godambe_matrix(
            params_0,
            model_func,
            model_args,
            means,
            varcovs,
            None,
            delta=delta,
            get_hessian=True
        )
        uncerts = np.sqrt(np.diag(np.linalg.inv(H)))

    elif method == "godambe":
        if bootstrap_reps is None:
            raise ValueError('we need bootstrap_reps to use `godambe` method!')
        G, _, __ = compute_godambe_matrix(
            params_0,
            model_func,
            model_args,
            means,
            varcovs,
            bootstrap_reps,
            delta=delta,
            get_hessian=False
        )
        uncerts = np.sqrt(np.diag(np.linalg.inv(G)))
    else:
        return

    return param_names, params_0, uncerts


def compute_godambe_matrix(
    params_0,
    model_func,
    model_args,
    means,
    varcovs,
    bootstrap_reps,
    delta=0.01,
    get_hessian=False,
    verbose=False
):
    """

    """
    def obj_func(params, means, varcovs, model_args):
        
        key = tuple(params)
        if key in _model_cache:
            model = _model_cache[key]
        else:
            model = model_func(params, model_args)
            _model_cache[key] = model

        return composite_ll(model, means, varcovs)

    H = -compute_hessian(
        params_0, 
        obj_func, 
        model_args,
        means,
        varcovs,
        delta=delta
    )
    if verbose:
        print(get_time(), "computed Hessian")

    if get_hessian:
        return H

    J = np.zeros((len(params_0), len(params_0)))
    for i, bootmeans in enumerate(bootstrap_reps):
        cU = compute_gradient(
            params_0, 
            obj_func, 
            model_args,
            bootmeans,
            varcovs,
            delta=delta
        )
        if verbose:
            print(get_time(), f"computed gradient for bootstrap set {i}")
        cJ = np.matmul(cU, cU.T)
        J += cJ
    J = J / len(bootstrap_reps)
    J_inv = np.linalg.inv(J)
    G = np.matmul(np.matmul(H,  J_inv), H)

    return G, H, J


def compute_hessian(p0, obj_func, model_args, means, varcovs, delta=0.01):
    """
    Compute the approximate Hessian matrix of the log-lik function. Uses 
    original data (not the bootstrap).
    """
    f0 = obj_func(p0, means, varcovs, model_args)
    hs = delta * p0
    H = np.zeros((len(p0), len(p0)), dtype=np.float64)
    for i in range(len(p0)):
        for j in range(i, len(p0)):
            p = np.array(p0, copy=True, dtype=np.float64)
            if i == j:
                p[i] = p0[i] + hs[i]
                fp = obj_func(p, means, varcovs, model_args)
                p[i] = p0[i] - hs[i]
                fm = obj_func(p, means, varcovs, model_args)
                element = (fp - 2 * f0 + fm) / hs[i] ** 2
            else:
                p[i] = p0[i] + hs[i]
                p[j] = p0[j] + hs[j]
                fpp = obj_func(p, means, varcovs, model_args)
                p[i] = p0[i] + hs[i]
                p[j] = p0[j] - hs[j]
                fpm = obj_func(p, means, varcovs, model_args)
                p[i] = p0[i] - hs[i]
                p[j] = p0[j] + hs[j]
                fmp = obj_func(p, means, varcovs, model_args)
                p[i] = p0[i] - hs[i]
                p[j] = p0[j] - hs[j]
                fmm = obj_func(p, means, varcovs, model_args)
                element = (fpp - fpm - fmp + fmm) / (4 * hs[i] * hs[j])
            H[i, j] = element
            H[j, i] = element

    return H


def compute_gradient(p0, obj_func, model_args, means, varcovs, delta=0.01):
    """
    """
    hs = delta * p0
    gradient = np.zeros((len(p0), 1))
    for i in range(len(p0)):
        p = np.array(p0, copy=True, dtype=float)
        p[i] = p0[i] + hs[i]
        fp = obj_func(p, means, varcovs, model_args)
        p[i] = p0[i] - hs[i]
        fm = obj_func(p, means, varcovs, model_args)
        gradient[i] = (fp - fm) / (2 * hs[i])

    return gradient


def indicator(n, i):
    # get an indicator vector of length n, with element i equal to 1
    if i >= n:
        raise ValueError("invalid index, must have i < n")
    arr = np.zeros(n, dtype=np.int64)
    arr[i] = 1

    return arr


## Utilities and printing functions


def graph_data_overlap(graph, pop_ids):
    """
    Find the populations which occur mutually in a Demes graph and a list of 
    population names.

    :param graph: Demes graph or the path and file name leading to a .yaml file
        specifying a demes graph. 
    """
    if isinstance(graph, str):
        graph = demes.load(graph)
    deme_names = [d.name for d in graph.demes]
    overlaps = [pop for pop in pop_ids if pop in deme_names]

    return overlaps


def get_time():
    """
    Get a string representing the date and time.
    """
    return "[" + datetime.strftime(datetime.now(), "%d-%m-%y %H:%M:%S") + "]"


## Caches


_model_cache = dict()






_out_of_bounds = -1e10
_n_calls = 0



def _moments_d_plus(
    graph,
    data=None,
    theta=None,
    rs=None,
    bins=None,
    u=None,
    approximation='trapezoid',
    sampled_demes=None,
    sample_times=None,
    phased=False
):

    methods = ['midpoint', 'trapezoid', 'Simpsons', None]
    if approximation not in methods: 
        raise ValueError(f"{approximation} is not a valid method")
    if u is None and theta is None:
        raise ValueError("you must provide `u` or `theta`")
    
    if isinstance(graph, str):
        graph = demes.load(graph)
    if data is not None:
        sampled_demes = data['pops']
        bins = data['bins']
    else:
        if sampled_demes is not None:
            graph_demes = [d.name for d in graph.demes]
            for d in sampled_demes:
                if d not in graph_demes: 
                    raise ValueError(f'deme {d} is not present in graph!')
        else:
            sampled_demes = [d.name for d in graph.demes if d.end_time == 0]
    if sample_times is None:
        end_times = {d.name: d.end_time for d in graph.demes}
        sample_times = [end_times[pop] for pop in sampled_demes]
    else:
        assert len(sample_times) == len(sampled_demes)
    if rs is None:
        if bins is None:
            raise ValueError("you must provide `bins` or `rs`")
        rs = get_rs(bins, approximation)   

    y = moments.Demes.LD(
        graph,
        sampled_demes,
        sample_times=sample_times,
        theta=theta,
        r=rs,
        u=u
    )
    num_demes = len(sampled_demes)
    indices = [(i, j) for i in range(num_demes) for j in range(i, num_demes)]
    raw_stats = np.zeros((len(rs), len(indices)))
    for k, (i, j) in enumerate(indices):
        if i == j:
            phasing = True
        else:
            phasing = phased
        raw_stats[:, k] = y.H2(i, j, phased=phasing)
    stats = approximate_Dplus(raw_stats, approximation)
    stats_H = np.vstack((stats, y.H()))
    model = {'means': stats_H, 'pops': sampled_demes, 'bins': bins}

    return model


def _get_rs(bins, approximation):

    key = (str(bins), approximation)
    if key in _rs_cache:
        rs = _rs_cache[key]

    elif approximation is None:
        rs = bins

    elif approximation == 'midpoint':
        rs = bins[:-1] + (bins[1:] - bins[:-1]) / 2

    elif approximation == 'trapezoid':
        rs = bins
    
    elif approximation == 'Simpsons':
        midpoints = (bins[1:] - bins[:-1]) / 2
        rs = np.sort(np.concatenate((bins, bins[:-1] + midpoints)))

    else:
        raise ValueError(f"{approximation} is not a valid method")

    return rs


def _approximate_Dplus(raw_stats, approximation):

    if approximation is None:
        ret = raw_stats

    elif approximation == 'midpoint':
        ret = raw_stats

    elif approximation == 'trapezoid':
        ret = 1/2 * (raw_stats[:-1] + raw_stats[1:])
    
    elif approximation == 'Simpsons':
        ret = (
            1/6 * raw_stats[:-1:2] 
            + 2/3 * raw_stats[1::2] 
            + 1/6 * raw_stats[2::2]
        )
    else:
        raise ValueError(f"{approximation} is not a valid method")

    return ret


def _load_statistics(file, graph=None):

    with open(file, 'rb') as fin:
        dic = pickle.load(fin)
    _data = dic[next(iter(dic))]
    if graph is not None:
        data = parsing.subset_statistics(_data, graph=graph)
    return data


## optimization functions


def __object_func(
    p,
    builder,
    options,
    data,
    u=None,
    lower_bounds=None,
    upper_bounds=None,
    constraints=None,
    verbose=None,
    one_locus=False
):
    
    global _n_calls
    _n_calls += 1
    if lower_bounds is not None and np.any(p < lower_bounds):
        return -_out_of_bounds
    elif upper_bounds is not None and np.any(p > upper_bounds):
        return -_out_of_bounds
    elif constraints is not None and np.any(constraints(p) <= 0):
        return -_out_of_bounds
    builder = Inference._update_builder(builder, options, p)
    graph = demes.Graph.fromdict(builder)
    model = moments_d_plus(graph, u=u, data=data)
    ll = compute_ll(model, data, one_locus=one_locus)
    if verbose > 0 and _n_calls % verbose == 0:
        print_status(_n_calls, ll, p)

    return -ll


def __object_func_log(logp, *args, **kwargs):
    
    p = np.exp(logp - 1)
    return _object_func(p, *args, **kwargs)


def _optimize(
    graph_file,
    param_file,
    data,
    u=None,
    method='fmin',
    max_iter=1000,
    max_calls=None,
    log=False,
    verbose=1,
    out_file=None,
    one_locus=False,
    perturb=False
):
    """
    Fit a graph defined in `graph_file` and parameterized by `param_file`
    to `data` using `objective_func` using a scipy optimization routine.
    """
    print(get_time(), f"fitting D+ to data for demes {data['pops']}")
    builder = Inference._get_demes_dict(graph_file)
    options = Inference._get_params_dict(param_file)
    pnames, p0, lower_bounds, upper_bounds = \
        Inference._set_up_params_and_bounds(options, builder)
    constraints = Inference._set_up_constraints(options, pnames)
    if u is None:
        raise ValueError("you must provide `u`")
    if perturb > 0: 
        p0 = Inference._perturb_params_constrained(
            p0, 
            perturb, 
            lower_bound=lower_bounds, 
            upper_bound=upper_bounds,
            cons=constraints
        )
        print_p0 = p0
    else:
        print_p0 = p0
    if log:
        objective = _object_func_log
        p0 = np.log(p0) + 1
    else:
        objective = _object_func
    print_status(0, 'pnames', pnames)
    print_status(0, 'p0', print_p0)
    
    warn = None
    args = (
        builder,
        options,
        data,
        u,
        lower_bounds,
        upper_bounds,
        constraints,
        verbose,
        one_locus
    )
    
    methods = ['fmin', 'powell', 'bfgs', 'lbfgsb']
    if method not in methods:
        raise ValueError(f"{method} is not a valid method")
    
    if method == 'fmin':
        output = scipy.optimize.fmin(
            objective,
            p0,
            args=args,
            maxiter=max_iter,
            maxfun=max_calls,
            full_output=True
        )
        popt, fopt, num_iter, func_calls, flag = output[:5]

    elif method == 'powell':
        output = scipy.optimize.fmin_powell(
            objective,
            p0,
            args=args,
            maxiter=max_iter,
            full_output=True,
        )
        popt, fopt, _, num_iter, func_calls, flag = output[:6]

    elif method == 'bfgs':
        if log:
            epsilon = 1e-3
        else:
            epsilon = None
        output = scipy.optimize.fmin_bfgs(
            objective,
            p0,
            args=args,
            maxiter=max_iter,
            epsilon=epsilon,
            disp=False,
            full_output=True
        )
        popt, fopt, _, __, func_calls, grad_calls, flag = output[:7]
        num_iter = grad_calls

    elif method == 'lbfgsb':
        if log:
            bounds = list(
                zip(np.log(lower_bounds) + 1, np.log(upper_bounds) + 1)
            )
            epsilon = 1e-5
            pgtol = 1e-5
        else:
            bounds = list(zip(lower_bounds, upper_bounds))
            epsilon = 1e-8
            pgtol = 1e-5
        output = scipy.optimize.fmin_l_bfgs_b(
            objective,
            p0,
            args=args,
            maxiter=max_iter,
            bounds=bounds,
            epsilon=epsilon,
            pgtol=pgtol,
            approx_grad=True
        )
        popt, fopt, output_dict = output
        num_iter = output_dict['nit']
        func_calls = output_dict['funcalls']
        flag = output_dict['warnflag']
        warn = output_dict["task"]

    else:
        return 1

    if log: 
        popt = np.exp(popt - 1)

    global _n_calls
    print_status(_n_calls, 'popt:', popt)
    info = dict(
        method=method,
        objective_func=objective.__name__,
        fopt=-fopt,
        max_iter=max_iter,
        num_iter=num_iter,
        func_calls=func_calls,
        flag=flag,
        warn=warn,
        u=u
    )
    print('\n'.join([f'\t{key}: {info[key]}' for key in info]))
    builder = Inference._update_builder(builder, options, popt)
    graph = demes.Graph.fromdict(builder)
    graph.metadata['opt_info'] = info

    if out_file is not None: 
        demes.dump(graph, out_file)
    else: 
        print(graph)
    
    return pnames, popt, fopt, graph


def _print_status(n_calls, ll, p):

    """
    Print the number of function calls, the log-likelihood, and the current 
    parameter values.
    """
    t = utils.get_time()
    _n = f'{n_calls:<4}'
    if isinstance(ll, float):
        _ll = f'{np.round(ll, 2):>10}'
    else:
        _ll = f'{ll:>10}'
    fmt_p = []
    for x in p:
        if isinstance(x, str):
            fmt_p.append(f'{x:>10}')
        else:
            if x > 1:
                fmt_p.append(f'{np.round(x, 1):>10}')
            else:
                sci = np.format_float_scientific(x, 2, trim='k')
                fmt_p.append(f'{sci:>10}')
    _p = ''.join(fmt_p)
    print(t, _n, _ll, '[', _p, ']')





