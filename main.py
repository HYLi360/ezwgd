from pprint import pprint
from ezwgd.evolution import dnds
from ezwgd.utils import parse


def main():
    genome = parse.Genome("test_dataset/testing.gff3", "test_dataset/testing.fna")

    print(genome.genelist)
if __name__ == "__main__":
    main()
