#!/usr/bin/env python3
# WHTOOLs MATCALIB 2026 - Prony to PRF (Power Law) Converter
"""
Convert Abaqus *VISCOELASTIC,TIME=PRONY parameters to PRF Power Law (FLAG_VISC=3).

Conversion formula:
    For each Prony term (g_i, k_i, tau_i) with instantaneous shear modulus mu:
        stiffn = g_i                               (stress fraction)
        A1     = 1 / (g_i * mu * tau_i)            (creep rate)
        EXPN   = 1                                 (linear stress dependence)
        EXPM   = 0                                 (no strain hardening)
        k_i    = 0                                 (bulk relaxation, typically zero)

Usage:
    python prony_to_prf.py --mu=1.0 0.5,1.0 0.3,0.1 --output ccx
    python prony_to_prf.py --c10=549.6 0.1986,2.81e-8 0.1828,2.81e-6 --solver optistruct
    python prony_to_prf.py --mu=1.0 --k=1000 0.5,1.0 --output both
"""
import argparse
import sys
from dataclasses import dataclass, field
from typing import List


@dataclass
class PronyTerm:
    g: float      # shear relaxation ratio
    tau: float    # relaxation time
    k: float = 0.0  # bulk relaxation ratio (typically 0)


@dataclass
class PRFTerm:
    stiffn: float
    flag_visc: int = 3     # Power Law
    A1: float = 0.0
    EXPN: float = 1.0       # linear stress
    EXPM: float = 0.0       # no strain hardening


@dataclass
class PRFConfig:
    n_network: int
    flag_he: int = 1        # Polynomial
    flag_pl: int = 0
    C10: float = 0.5
    mu: float = 1.0          # instantaneous shear modulus
    K: float = 1000.0        # bulk modulus (CCX convention: K=2*D1)
    terms: List[PRFTerm] = field(default_factory=list)


def prony_to_prf(prony_terms: List[PronyTerm], mu: float, K: float = 1000.0,
                 C10: float = None) -> PRFConfig:
    """Convert Prony series terms to PRF Power Law configuration."""
    if C10 is None:
        C10 = mu / 2.0
    n = len(prony_terms)
    terms = []
    for pt in prony_terms:
        if pt.g < 1e-20:
            continue
        A1 = 1.0 / (pt.g * mu * pt.tau)
        terms.append(PRFTerm(stiffn=pt.g, A1=A1))
    return PRFConfig(n_network=n, C10=C10, mu=mu, K=K, terms=terms)


def build_ccx_constants(cfg: PRFConfig) -> list:
    """Generate CCX *USER MATERIAL constants array."""
    # D1 = K/2 (CCX convention)
    D1 = cfg.K / 2.0
    c = [
        float(cfg.n_network), float(cfg.flag_he), float(cfg.flag_pl),
        cfg.C10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        D1, 0.0, 0.0,
        1.0, cfg.K,    # IHYPER=1, RBULK_BACKUP
    ]
    for t in cfg.terms:
        c += [t.stiffn, float(t.flag_visc), t.A1, t.EXPN, t.EXPM]
    c += [1.0, cfg.K]  # G, RBULK_FINAL
    return c


def format_8perline(constants: list) -> str:
    """Format constants exactly 8 per line, last line padded with 0.0."""
    all_vals = list(constants) + [0.0]
    lines = []
    for i in range(0, len(all_vals), 8):
        chunk = all_vals[i:i+8]
        if not chunk:
            continue
        while len(chunk) < 8:
            chunk.append(0.0)
        lines.append(', '.join(f'{v:.15g}' for v in chunk))
    return '\n'.join(lines)


def calc_depvar(cfg: PRFConfig) -> int:
    return 2 + 12 * cfg.n_network


def main():
    parser = argparse.ArgumentParser(
        description='Convert Prony series to PRF Power Law parameters')
    parser.add_argument('terms', nargs='+', help='Prony terms: g,tau or g,k,tau')
    parser.add_argument('--mu', type=float, default=None,
                       help='Instantaneous shear modulus (mu = 2*C10)')
    parser.add_argument('--c10', type=float, default=None,
                       help='Neo-Hookean C10 (alternative to --mu)')
    parser.add_argument('--k', type=float, default=1000.0,
                       help='Bulk modulus (CCX convention: K=2*D1, default: 1000)')
    parser.add_argument('--output', choices=['ccx', 'optistruct', 'both'], default='both',
                       help='Output format (default: both)')
    parser.add_argument('--flag-he', type=int, default=1,
                       help='FLAG_HE type (default: 1=Polynomial)')
    parser.add_argument('--noproty', type=int, default=None,
                       help='Override N_NETWORK (default: number of Prony terms)')
    args = parser.parse_args()

    # Parse Prony terms
    prony_terms = []
    for t_str in args.terms:
        parts = t_str.split(',')
        if len(parts) == 2:
            g, tau = float(parts[0]), float(parts[1])
            prony_terms.append(PronyTerm(g=g, tau=tau))
        elif len(parts) == 3:
            g, k, tau = float(parts[0]), float(parts[1]), float(parts[2])
            prony_terms.append(PronyTerm(g=g, k=k, tau=tau))
        else:
            print(f"ERROR: Invalid term '{t_str}'. Use 'g,tau' or 'g,k,tau'")
            sys.exit(1)

    if not prony_terms:
        print("ERROR: No valid Prony terms provided")
        sys.exit(1)

    # Determine mu
    if args.mu is not None:
        mu = args.mu
        c10 = mu / 2.0
    elif args.c10 is not None:
        c10 = args.c10
        mu = 2.0 * c10
    else:
        print("ERROR: Specify --mu or --c10")
        sys.exit(1)

    # Convert
    cfg = prony_to_prf(prony_terms, mu, args.k, c10)
    if args.noproty is not None:
        cfg.n_network = args.noproty
    cfg.flag_he = args.flag_he

    constants = build_ccx_constants(cfg)
    depvar = calc_depvar(cfg)
    sum_stiffn = sum(t.stiffn for t in cfg.terms)
    eq_frac = 1.0 - sum_stiffn

    # Print conversion summary
    print("=" * 70)
    print("  Prony -> PRF Power Law Conversion")
    print("=" * 70)
    print(f"  mu = {mu:.6g}  (C10 = {c10:.6g})")
    print(f"  K  = {cfg.K:.6g}  (D1 = {cfg.K/2:.6g})")
    print(f"  Equilibrium fraction: {eq_frac:.6f}  (1 - sum(g_i))")
    print()

    for i, (pt, prf_t) in enumerate(zip(prony_terms, cfg.terms)):
        print(f"  Prony term {i+1}: g={pt.g:.6g}, k={pt.k:.6g}, tau={pt.tau:.6g}")
        print(f"  -> PRF Power Law: stiffn={prf_t.stiffn:.6g}, "
              f"A1={prf_t.A1:.6g}, EXPN={prf_t.EXPN}, EXPM={prf_t.EXPM}")
    print()

    # CCX output
    if args.output in ('ccx', 'both'):
        print("=" * 70)
        print("  CalculiX *USER MATERIAL")
        print("=" * 70)
        print(f"  N_NETWORK = {cfg.n_network}")
        print(f"  DEPVAR    = {depvar}")
        print(f"  CONSTANTS = {len(constants)}")
        print()
        print(f"*USER MATERIAL,CONSTANTS={len(constants)}")
        print(format_8perline(constants))
        print(f"*DEPVAR")
        print(f"{depvar}")
        print()

    # OptiStruct output
    if args.output in ('optistruct', 'both'):
        print("=" * 70)
        print("  OptiStruct MATUSR")
        print("=" * 70)
        print(f"  NDEPVAR   = {depvar}")
        print(f"  CONSTANTS = {len(constants)}")
        print()
        const_str = ', '.join(f'{v:.15g}' for v in constants)
        print(f"MATUSR, <MID>, <USUBID>, {depvar}, <GROUP>, <DENSITY>, PROPERTY,")
        # Split long line
        for i in range(0, len(const_str), 120):
            print(f"        {const_str[i:i+120]}")
        print()


if __name__ == '__main__':
    main()
