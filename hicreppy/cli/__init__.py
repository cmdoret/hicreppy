# -*- coding: utf-8 -*-
import sys
import click
import cooler as clr
from .. import hicrep as cr
from .. import __version__


@click.command()
@click.option(
    "--h-max",
    "-r",
    default=10,
    show_default=True,
    help="Maximum value of the smoothing parameter h to explore. All consecutive integer values from 0 to this value will be tested.",
)
@click.option(
    "--max-dist",
    "-m",
    default=100000,
    show_default=True,
    help="Maximum distance at which to compute the SCC, in basepairs.",
)
@click.option(
    "--blacklist",
    "-b",
    default="",
    help="Exclude those chromosomes in the analysis. List of comma-separated chromosome names.",
)
@click.option(
    "--whitelist",
    "-w",
    default="",
    help="Only include those chromosomes in the analysis. List of comma-separated chromosome names.",
)
@click.argument("cool1", type=click.Path(exists=False))
@click.argument("cool2", type=click.Path(exists=False))
def htrain_cmd(cool1, cool2, max_dist, h_max, whitelist, blacklist):
    """Find the optimal value for smoothing parameter h.
    The optimal h-value is printed to stdout. Run informations are printed to stderr.
    """
    c1 = clr.Cooler(cool1)
    c2 = clr.Cooler(cool2)
    wl = _parse_cli_chroms(whitelist)
    bl = _parse_cli_chroms(blacklist)
    # Print best h value to stdout
    print(cr.h_train(c1, c2, max_dist, h_max, whitelist=wl, blacklist=bl))


@click.command()
@click.option(
    "--h-value",
    "-v",
    default=10,
    show_default=True,
    help="Value of the smoothing parameter h to use. Should be an integer value >= 0.",
)
@click.option(
    "--max-dist",
    "-m",
    default=100000,
    show_default=True,
    help="Maximum distance at which to compute the SCC, in basepairs.",
)
@click.option(
    "--subsample",
    "-s",
    default=None,
    show_default=True,
    help="Subsample contacts from both matrices to target value. Leave to 0 to disable subsampling.",
)
@click.option(
    "--blacklist",
    "-b",
    default="",
    help="Exclude those chromosomes in the analysis. List of comma-separated chromosome names.",
)
@click.option(
    "--whitelist",
    "-w",
    default="",
    help="Only include those chromosomes in the analysis. List of comma-separated chromosome names.",
)
@click.argument("cool1", type=click.Path(exists=False))
@click.argument("cool2", type=click.Path(exists=False))
def genome_scc_cmd(
    cool1, cool2, max_dist, h_value, subsample, whitelist, blacklist
):
    """Compute the stratum-adjusted correlation coefficient for input matrices"""
    c1 = clr.Cooler(cool1)
    c2 = clr.Cooler(cool2)
    wl = _parse_cli_chroms(whitelist)
    bl = _parse_cli_chroms(blacklist)

    # Make sure a proper value was given to subsample
    if subsample is not None:
        try:
            subsample = float(subsample)
        except ValueError:
            raise(
                "Subsample must be a float between 0 and 1 or an integer."
            )
    print(
        cr.genome_scc(
            c1,
            c2,
            max_dist,
            h_value,
            subsample=subsample,
            whitelist=wl,
            blacklist=bl,
        )
    )


@click.group()
@click.version_option(version=__version__)
def cli():
    """
    Type -h or --help after any subcommand for more information.
    """
    ...


cli.add_command(htrain_cmd, name="htrain")
cli.add_command(genome_scc_cmd, name="scc")


def _parse_cli_chroms(chroms):
    """Process a string of comma separated chromosomes into a list"""
    if chroms == "":
        return None
    else:
        return chroms.split(",")
