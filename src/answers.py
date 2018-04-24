from sortedcontainers import SortedList
from structs import Segment, Pixel, SweepLine as SL, Point, normalize, get_pixel, add_pixel_to_seg
from common import hot, current, line, segments



def bcheck(point):
	pixel = Pixel(point)
	return (pixel.is_on_top(point) and pixel.get_top_neighbour(point) in current) or (pixel.is_on_bottom(point) and pixel.get_bottom_neighbour(point) in current)


def segment_endpoint_answer(point, segment):
	pixel = get_pixel(point)
	if pixel not in current:
		# до этого в пикселе не было обнаружено критических точек
		heat_answer(pixel)
	# если событие -- начало отрезка
	if point == segment.start:
		add_pixel_to_seg(pixel, segment)
		# игнорируем ту часть отрезка, что лежит внутри пикселя
		for intersection in segment.intersections(pixel):
			if intersection is not None:
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break;
	# если событие -- конец отрезка
	else:
		line.remove(segment)

def segseg_intersection_answer(point):
	new_point = normalize(point)
	pixel = get_pixel(new_point)
	if bcheck(point):
		return
	if pixel not in current:
		heat_answer(pixel)

def segpix_intersection_answer(point, segment, pixel):
	# добавляем отрезок в соответствующий список
	add_pixel_to_seg(pixel, segment)				
	# игнорируем фрагмент отрезка, находящийся внутри пикселя
	line.remove(segment)
	lastinter = None
	for intersection in segment.intersections(pixel):
		if ((intersection is not None) and (intersection > point)) and ((lastinter is None) or (intersection > lastinter)):
			lastinter = intersection
	if (lastinter is not None):
		line.push(SL.Event(SL.Event.Type.SEG_REINSERT, lastinter, segment=segment, pixel=pixel))


def segment_reinsertion_answer(point, segment, pixel):
	line.insert(segment)
	neighbour = pixel.get_neighbour(point)
	if (neighbour is not None) and (neighbour in current):
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
	if (pixel in current):
		return
	current.add(pixel)
	l = []
	for s in segments:
		for intersection in s.intersections(pixel) :
			if (intersection is not None):
				add_pixel_to_seg(pixel, s)	
				if (intersection.x / intersection.z < line.xpos):
					l.append(s)
					break

	for segment in l:
		for intersection in segment.intersections(pixel):
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				line.remove(segment)
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break	
	top, bottom = pixel.top(), pixel.bottom()
	line.insert(top)
	line.push(SL.Event(SL.Event.Type.PIX_END, top.end, pixel=pixel))
	line.insert(bottom)
	line.push(SL.Event(SL.Event.Type.PIX_END, bottom.end, pixel=pixel))
