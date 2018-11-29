#!/usr/local/bin/python

from pyliftover import LiftOver
import os, sys

lo = LiftOver(os.getcwd() + '/' + sys.argv[1])
#print(sys.argv[1], sys.argv[2])
coord = lo.convert_coordinate(str(sys.argv[2]), int(sys.argv[3]))
print(coord)
#print(coord[0][0], coord[0][1])
