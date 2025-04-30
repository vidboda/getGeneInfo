import sys
from liftover import get_lifter

converter = get_lifter(str(sys.argv[1]), str(sys.argv[2]), one_based=True)
chrom = str(sys.argv[3])
pos = int(sys.argv[4])
print(converter[chrom][pos])
