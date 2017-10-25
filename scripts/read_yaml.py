#! /usr/bin/env python

import sys
import yaml

def main(argv):
    with open(argv[0], 'r') as fin:
        try:
            print(yaml.load(fin))
            return 0

        except yaml.YAMLError as err:
            print(err)
            return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
