from GH import Segment

s1 = Segment(Point(0, 0), Point(4, 4))
s2 = Segment(Point(0, 2), Point(4, 2))
print(s1.intersects(s2))