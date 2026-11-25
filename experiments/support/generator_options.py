"""Composition-search CLI options shared by the size and runtime experiments."""
from itertools import product

COUNT_DEFAULTS = {
    "composition-orbit-rounds": 12,
    "composition-orbit-min-gain": 1,
    "composition-orbit-patience": 3,
    "composition-orbit-random-seed": 923407,
    "composition-orbit-pair-visits": 0,
    "composition-orbit-seed-terms": 65536,
    "composition-orbit-product-terms": 65536,
    "composition-nilpotent-rounds": 8,
    "composition-certificate-visits": 8192,
    "composition-nilpotent-seeds": 8,
    "composition-nilpotent-products": 8192,
    "composition-nilpotent-pair-visits": 131072,
    "composition-nilpotent-seed-basis": 8,
    "composition-nilpotent-seed-terms": 256,
    "composition-nilpotent-product-terms": 65536,
}


def nonnegative_count(value):
    if not value.isascii() or not value.isdecimal() or int(value) > 2**64 - 1:
        raise ValueError("Expected a nonnegative 64-bit integer")
    return int(value)


def add_generator_arguments(parser, *, rounds=None, modes=None, defaults=None):
    for name, default in COUNT_DEFAULTS.items():
        default = (defaults or {}).get(name, default)
        settings = dict(type=nonnegative_count, default=default,
                        help=f"Composition search count (default: {default})")
        if name == "composition-orbit-rounds" and rounds is not None:
            settings.update(nargs="+", default=rounds, help="Orbit round limits to compare (default: 3 6 12)")
        parser.add_argument("--" + name, **settings)
    settings = dict(choices=("sum", "random", "hybrid"), default="hybrid",
                    help="Orbit seeds: all remaining maps, random half-subsets, or alternating broad/random")
    if modes is not None:
        settings.update(nargs="+", default=modes)
    parser.add_argument("--composition-orbit-seed-mode", **settings)


def generator_settings(args):
    return {name: getattr(args, name.replace("-", "_"))
            for name in (*COUNT_DEFAULTS, "composition-orbit-seed-mode")}


def generator_sweep(args):
    settings = generator_settings(args)
    rounds = settings["composition-orbit-rounds"]
    modes = settings["composition-orbit-seed-mode"]
    return {f"generators-{mode}-{count}": dict(settings, **{
                "composition-orbit-rounds": count, "composition-orbit-seed-mode": mode})
            for mode, count in product(dict.fromkeys(modes), dict.fromkeys(rounds))}


def generator_flags(settings):
    return [f"--{name}={value}" for name, value in settings.items()]
