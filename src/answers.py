from sortedcontainers import SortedList
from structs import Segment, Pixel, SweepLine as SL, Point, normalize
from common import hot, current, line, segments

heated = set()

def scheck(point):
	new_point = normalize(point)
	if ((new_point.x, new_point.y, new_point.z) in heated):
		return True
	else:
		return False

def bcheck(point):
	pixel = Pixel(point)
	if pixel.is_on_top(point):
		print("ontop")
		print(pixel.get_top_neighbour(point).center)
		if scheck(pixel.get_top_neighbour(point).center):
			return True
	if pixel.is_on_bottom(point):
		print("onbottom")
		if scheck(pixel.get_bottom_neighbour(point).center):
			return True
	if pixel.is_on_left(point):
		print("onleft")
		if scheck(pixel.get_left_neighbour(point).center):
			return True
	if pixel.is_on_right(point):
		print("onright")
		if scheck(pixel.get_right_neighbour(point).center):
			return True
	return False


def segment_endpoint_answer(point, segment):
	pixel = Pixel(point)
	print("rounded ", pixel.center)
	if pixel not in current:
		print("nice")
		# до этого в пикселе не было обнаружено критических точек
		heat_answer(pixel)
	# если событие -- начало отрезка
	if point == segment.start:
		pixel.segs.append(segment)
		# игнорируем ту часть отрезка, что лежит внутри пикселя
		intersection = segment.intersects(pixel)
		if intersection is not None:
			line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
	# если событие -- конец отрезка
	else:
		line.remove(segment)

def segseg_intersection_answer(point):
	new_point = normalize(point)
	pixel = Pixel(new_point)
	if bcheck(point):
		print(new_point, " hurray	")
		return
	if pixel not in current:
		print("from segseg ", new_point)
		for a in heated:
			print(a)
		print("end")
		heat_answer(pixel)

def segpix_intersection_answer(point, segment, pixel):
	# добавляем отрезок в соответствующий список
	if pixel.is_on_top(point):
		pixel.upper.add(segment)
	else:
		pixel.lower.add(segment)
	pixel.segs.append(segment)
	# игнорируем фрагмент отрезка, находящийся внутри пикселя
	line.remove(segment)
	# TODO FIXME INTERSECTS TO THE RIGHT OF POINT
	intersection = segment.intersects(pixel)
	if intersection is not None and (intersection.x/intersection.z) > (point.x/point.z):
		line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))

def segment_reinsertion_answer(point, segment, pixel):
	print("inserting ", segment)
	line.insert(segment)
	if pixel.is_on_top(point):
		pixel.upper.add(segment)
	elif pixel.is_on_bottom(point):
		pixel.lower.add(segment)

	neighbour = pixel.get_neighbour(point)
	if neighbour in current:
		segpix_intersection_answer(point, segment, neighbour)

def pixel_end_answer(point, pixel):
	if pixel.is_on_top(point):
		line.remove(pixel.top())
	else:
		line.remove(pixel.bottom())
	hot.extend(current)
	current.clear()

def heat_answer(pixel):
	pixel.center = normalize(pixel.center)
	if ((pixel.x, pixel.y, pixel.z) in heated):
		return
	heated.add((pixel.x, pixel.y, pixel.z))
	current.add(pixel)
	print("hhhhhhhheating ", pixel.center)
	l = []
	for s in segments:
		intersections = [s.intersects(pixel.top()), s.intersects(pixel.bottom()), s.intersects(pixel.right()), s.intersects(pixel.left())]
		for intersection in intersections :
			if (intersection is not None and (intersection.x/intersection.z) < line.xpos):
				l.append(s)
				break

	for segment in l:
		pixel.segs.append(segment)
		intersections = [segment.intersects(pixel.top()), segment.intersects(pixel.bottom()), segment.intersects(pixel.right()), s.intersects(pixel.left())]
		for intersection in intersections:
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				line.remove(segment)
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break	
	# pix.build()
	print ("inserting pixel")
	top, bottom = pixel.top(), pixel.bottom()
	print (top.end)
	print (bottom.end)
	print ("____________")
	line.insert(top)
	line.push(SL.Event(SL.Event.Type.PIX_END, top.end, pixel=pixel))
	line.insert(bottom)
	line.push(SL.Event(SL.Event.Type.PIX_END, bottom.end, pixel=pixel))
