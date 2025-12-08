"""
Cli of python application 
"""

from argparse import ArgumentParser
from typing import List 


def create_parser() -> ArgumentParser:
    """Create parser"""
    parser = ArgumentParser(
        prog = "Python Template",
        description= "Example python project"
    )
    # Add arguments here such as:
    #   parser.add_argument('posarg', choices=["foo", "bar"])
    #   parser.add_argument('-v', '--value', type=int)
    #   parser.add_argument('-f', '--flag', dest="myflag", action="store_true")
    #   parser.add_argument('-nf', '--no-flag', dest="myflag", action="store_false")
    return parser

def main(args: List[str] = None):
    parser = create_parser()
    args = parser.parser_args(args)

    # Do whatever the application does, such as:
    #   myfunc(myarg=args.value, is_flag=args.myflag)
    ...