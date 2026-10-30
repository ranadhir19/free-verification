"""The Section 9 bar (PREREGISTRATION_REBUILD.md s.5.6).

Imported operating point p_c = 0.0882; schedule p_c(r) = 0.0882 + kappa*(1-r).
Static bar: (1-rho*) = odds(p_c)/odds(p).
Model-conditional bar: solve odds(p)*(1-rho) = odds(p_c(r)) with
r = (1-p) + (1-rho)*p (sound constraint, sigma = 1), bisection to 1e-12.
"""
import math


def odds(p):
    return p / (1.0 - p)


def static_bar(p, pc=0.0882):
    omr = odds(pc) / odds(p)          # 1 - rho*
    rho = 1.0 - omr
    return {"rho_star": rho, "bits": math.log2(1.0 / omr)}


def model_bar(p, kappa=0.038, pc0=0.0882):
    # f(rho) = odds(p)*(1-rho) - odds(pc0 + kappa*(1 - r(rho))), r(rho)=(1-p)+(1-rho)*p
    def f(rho):
        omr = 1.0 - rho
        r = (1.0 - p) + omr * p
        return odds(p) * omr - odds(pc0 + kappa * (1.0 - r))
    lo, hi = 0.0, 1.0 - 1e-15
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-13:
            break
    rho = 0.5 * (lo + hi)
    omr = 1.0 - rho
    r = (1.0 - p) + omr * p
    pc = pc0 + kappa * (1.0 - r)
    return {"rho_star": rho, "r": r, "pc_prime": pc,
            "bits": math.log2(1.0 / omr)}


if __name__ == "__main__":
    # Self-check against the PRINTED record (old counts 3070/15564).
    p_old = 3070 / 15564
    print(f"printed p = {p_old:.6f} (printed 0.1973)")
    s = static_bar(p_old)
    print(f"static: rho*={s['rho_star']:.4f} bits={s['bits']:.3f}  (printed 0.6065 / 1.345)")
    for k in (0.0, 0.019, 0.038, 0.057):
        m = model_bar(p_old, k)
        print(f"kappa={k:<6} rho*={m['rho_star']:.4f} r={m['r']:.4f} "
              f"pc'={m['pc_prime']:.4f} bits={m['bits']:.3f}")
    print("printed: k=0 -> 1.345; k=0.019 -> 1.306; k=0.038 -> rho* 0.5845, "
          "r 0.885, pc' 0.0926, 1.267; k=0.057 -> 1.233")
