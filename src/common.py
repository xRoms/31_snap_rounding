from sortedcontainers import SortedList
from utils import bound
from random import uniform
from structs import rounddd, Point, Pixel, Segment, SweepLine as SL, mkpoint
import numpy as np

source = [ [(0.3, 0.3), (2.1, 4.4)], [(0.6, 2.0), (2.7, 2.)] ]
segments2  = list(map(lambda s: Segment(
		mkpoint(s[0][0], s[0][1]),
		mkpoint(s[1][0], s[1][1])
	), source))


segments = [Segment(mkpoint(4.1, 4.1), mkpoint(6.1, 0.1)), Segment(mkpoint(2.1, 2.1), mkpoint(8.1, 1.1)), 
			Segment(mkpoint(2, 0), mkpoint(6.1, 4)), Segment(mkpoint(2.1, -0.1), mkpoint(6.0, 3.9)),
			Segment(mkpoint(2, 1), mkpoint(3, 0))]

c = mkpoint(2.3, 4.5)
print(rounddd(c))

def generate_segs(n):
	rand = lambda : uniform(0, bound)
	build_segs([ [(rand(), rand()), (rand(), rand())] for i in range(n) ])

def build_segs(source):
	global segments
	segments.extend(list(map(lambda s: Segment(
		mkpoint(s[0][0], s[0][1]), 
		mkpoint(s[1][0], s[1][1])
	), source)))

hot = []
current = SortedList()
line = SL(segments)
