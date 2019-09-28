package datastruct;

import vecmath.Point3d;

/* Simple class for forming and manipulating kD-trees
 */
public class KDNode {

	// The amount of space that needs to be empty before a
	// "vacuum cut" (one half completely empty) is
	// considered
	private static final double VACUUM_RATIO = 0.5;

	// Maximum number of points in a leaf node
	private static final int MAX_LEAVES = 8;
	private final Point3d[] points;
	private int dimension; // 0=x, 1=y, 2=z, (-1=no cut)
	private double cutoff;
	private KDNode firstChild;
	private KDNode secondChild;
	private int firstLeaf;
	private int lastLeaf;

	// Construktor with points, point indices and corner coordinates
	public KDNode(Point3d[] points) {

		// System.out.println("tree initialized!");
		dimension = -1;
		cutoff = 0.0;
		firstChild = null;
		secondChild = null;
		this.points = points;
		firstLeaf = -1; // no points
		lastLeaf = -1;
	}

	private static boolean isSortedZ(Point3d[] points, int left, int right) {
		for (int i = left + 1; i <= right; i++) {
			if (points[i].z < points[i - 1].z)
				return false;
		}
		return true;
	}

	public void build(int first, int last, double cx0, double cx1, double cy0, double cy1, double cz0, double cz1) {

		double x0, x1, y0, y1, z0, z1;

		// System.out.println("Build with nodes "+first+" to "+last);

		// build tree
		if ((last - first + 1) <= MAX_LEAVES) {
			firstLeaf = first;
			lastLeaf = last;
			// System.out.println(" -> leaf node");
			// finished
		} else {
			// check "vacuum" dimensions
			// current bounding box
			x0 = cx1;
			x1 = cx0;
			y0 = cy1;
			y1 = cy0;
			z0 = cz1;
			z1 = cz0;
			for (int i = first; i <= last; i++) {
				Point3d p = points[i];
				if (p.x < x0)
					x0 = p.x;
				if (p.x > x1)
					x1 = p.x;
				if (p.y < y0)
					y0 = p.y;
				if (p.y > y1)
					y1 = p.y;
				if (p.z < z0)
					z0 = p.z;
				if (p.z > z1)
					z1 = p.z;
			}
			// assert that all points are within our bounding box!
			if (cx0 > x0) {
				System.err.println("X min out of bounds!");
			}
			if (cx1 < x1) {
				System.err.println("X max out of bounds!");
			}
			if (cy0 > y0) {
				System.err.println("Y min out of bounds!");
			}
			if (cy1 < y1) {
				System.err.println("Y max out of bounds!");
			}
			if (cz0 > z0) {
				System.err.println("Z min out of bounds!");
			}
			if (cz1 < z1) {
				System.err.println("Z max out of bounds!");
			}
			// two cutoff-checks per dimension
			int which;
			double val;
			which = 0;
			val = (x0 - cx0) / (cx1 - cx0);
			if ((cx1 - x1) / (cx1 - cx0) > val) {
				which = 1;
				val = (cx1 - x1) / (cx1 - cx0);
			}
			if ((y0 - cy0) / (cy1 - cy0) > val) {
				which = 2;
				val = (y0 - cy0) / (cy1 - cy0);
			}
			if ((cy1 - y1) / (cy1 - cy0) > val) {
				which = 3;
				val = (cy1 - y1) / (cy1 - cy0);
			}
			if ((z0 - cz0) / (cz1 - cz0) > val) {
				which = 4;
				val = (z0 - cz0) / (cz1 - cz0);
			}
			if ((cz1 - z1) / (cz1 - cz0) > val) {
				which = 5;
				val = (cz1 - z1) / (cz1 - cz0);
			}
			if (val > VACUUM_RATIO) {
				// vacuum cutoff!

				switch (which) {
					case 0 : // x min cutoff
						// System.out.println(" -> x min vacuum");
						dimension = 0; // x
						cutoff = x0;
						secondChild = new KDNode(points);
						secondChild.build(first, last, x0, cx1, cy0, cy1, cz0, cz1);
						firstChild = null;
						break;
					case 1 : // x max cutoff
						// System.out.println(" -> x max vacuum");
						dimension = 0; // x
						cutoff = x1;
						secondChild = null;
						firstChild = new KDNode(points);
						firstChild.build(first, last, cx0, x1, cy0, cy1, cz0, cz1);
						break;
					case 2 : // y min cutoff
						// System.out.println(" -> y min vacuum");
						dimension = 1; // y
						cutoff = y0;
						secondChild = new KDNode(points);
						secondChild.build(first, last, cx0, cx1, y0, cy1, cz0, cz1);
						firstChild = null;
						break;
					case 3 : // y max cutoff
						// System.out.println(" -> y max vacuum");
						dimension = 1; // y
						cutoff = y1;
						secondChild = null;
						firstChild = new KDNode(points);
						firstChild.build(first, last, cx0, cx1, cy0, y1, cz0, cz1);
						break;
					case 4 : // z min cutoff
						// System.out.println(" -> z min vacuum");
						dimension = 2; // z
						cutoff = z0;
						secondChild = new KDNode(points);
						secondChild.build(first, last, cx0, cx1, cy0, cy1, z0, cz1);
						firstChild = null;
						break;
					case 5 : // z max cutoff
						// System.out.println(" -> z max vacuum");
						dimension = 2; // z
						cutoff = z1;
						secondChild = null;
						firstChild = new KDNode(points);
						firstChild.build(first, last, cx0, cx1, cy0, cy1, cz0, z1);
						break;
				}

			} else {
				// no vacuum cutoff

				// get axis of greatest variance
				x0 = 0.0;
				x1 = 0.0;
				y0 = 0.0;
				y1 = 0.0;
				z0 = 0.0;
				z1 = 0.0;

				for (int i = first; i <= last; i++) {
					Point3d p = points[i];
					x0 += p.x;
					x1 += p.x * p.x;
					y0 += p.y;
					y1 += p.y * p.y;
					z0 += p.z;
					z1 += p.z * p.z;
				}
				x1 = x1 * (last - first + 1) - x0 * x0;
				y1 = y1 * (last - first + 1) - y0 * y0;
				z1 = z1 * (last - first + 1) - z0 * z0;

				int break_node;

				if (x1 > y1) {
					// x axis
					// System.out.println(" -> x axis split");
					// KHK
					// KHK war der Originalaufruf: sort(first, last, 0);
					sortKH(first, last, 0);
					break_node = (first + last) / 2;
					dimension = 0;
					cutoff = points[break_node].x;
					firstChild = new KDNode(points);
					firstChild.build(first, break_node, cx0, cutoff, cy0, cy1, cz0, cz1);
					secondChild = new KDNode(points);
					secondChild.build(break_node + 1, last, cutoff, cx1, cy0, cy1, cz0, cz1);
				} else if (y1 > z1) {
					// y axis
					// System.out.println(" -> y axis split");
					// KHK
					// KHK war der Originalaufruf: sort(first, last, 1);
					sortKH(first, last, 1);
					break_node = (first + last) / 2;
					dimension = 1;
					cutoff = points[break_node].y;
					firstChild = new KDNode(points);
					firstChild.build(first, break_node, cx0, cx1, cy0, cutoff, cz0, cz1);
					secondChild = new KDNode(points);
					secondChild.build(break_node + 1, last, cx0, cx1, cutoff, cy1, cz0, cz1);
				} else {
					// z axis
					// System.out.println(" -> z axis split");
					// KHK
					// KHK war der Originalaufruf: sort(first, last, 2);
					sortKH(first, last, 2);
					break_node = (first + last) / 2;
					dimension = 2;
					cutoff = points[break_node].z;
					firstChild = new KDNode(points);
					firstChild.build(first, break_node, cx0, cx1, cy0, cy1, cz0, cutoff);
					secondChild = new KDNode(points);
					secondChild.build(break_node + 1, last, cx0, cx1, cy0, cy1, cutoff, cz1);
				}

			}
		}

	}

	/*
	 * KHK: I had to replace the classic Quicksort implementation with this
	 * variation because of stack overflow errors.
	 * 
	 * 
	 * I got the idea for my interpretation of Quicksort from here:
	 * 
	 * // from: http://algs4.cs.princeton.edu/23quicksort/Quick3way.java.html //
	 * quicksort the arraypoints] using 3-way partitioning
	 * 
	 * public static void sort(Point3d[] points, int dim) { sort(points, 0,
	 * points.length - 1); switch(dim) { case 0: assert isSortedX(points); break;
	 * case 1: assert isSortedY(points); break; case 2: assert isSortedZ(points); }
	 * }
	 * 
	 * // quicksort the subarraypointsleft .. right] using 3-way partitioning
	 * private static void sort(Point3d[] points, int left, int right, int dim) {
	 */

	public int findNearest(Point3d point, double distance, double dx, double dy, double dz) {
		double d;
		int p = -1, p2 = -1;
		// all distances are squared!
		if (firstLeaf >= 0) {
			// Node is leaf node
			for (int i = firstLeaf; i <= lastLeaf; i++) {
				d = point.distanceSquared(points[i]);
				if (d < distance) {
					p = i;
					distance = d;
				}
			}
		} else {
			// composite node
			// // complete search for now!
			// if (firstChild != null) {
			// p = firstChild.findNearest(point, distance, cube_distance);
			// if (p >= 0) distance = point.distanceSquared(points[p]);
			// }
			// if (secondChild != null) {
			// int p2 = secondChild.findNearest(point, distance, cube_distance);
			// if (p2 >= 0) {
			// distance = point.distanceSquared(points[p2]);
			// p = p2;
			// }
			// }
			double splitSide = 0.0;
			switch (dimension) {
				case 0 :
					splitSide = cutoff - point.x;
					break;
				case 1 :
					splitSide = cutoff - point.y;
					break;
				case 2 :
					splitSide = cutoff - point.z;
					break;
			}
			// If splitSide is >= 0, the first child is nearer
			if (splitSide >= 0.0) {
				if (firstChild != null) {
					p = firstChild.findNearest(point, distance, dx, dy, dz);
				}
				if (p >= 0) {
					// success!
					// diminish search radius
					distance = point.distanceSquared(points[p]);
				}
				// get new cube distance
				switch (dimension) {
					case 0 :
						dx = (cutoff - point.x) * (cutoff - point.x);
						break;
					case 1 :
						dy = (cutoff - point.y) * (cutoff - point.y);
						break;
					case 2 :
						dz = (cutoff - point.z) * (cutoff - point.z);
						break;
				}
				d = dx + dy + dz;
				if ((d < distance) && (secondChild != null)) {
					// look only if child is closer than radius
					p2 = secondChild.findNearest(point, distance, dx, dy, dz);
				}
				if (p2 >= 0)
					p = p2;
			} else {
				// second child is nearer
				if (secondChild != null) {
					p = secondChild.findNearest(point, distance, dx, dy, dz);
				}
				if (p >= 0) {
					// success!
					// diminish search radius
					distance = point.distanceSquared(points[p]);
				}
				// get new cube distance
				switch (dimension) {
					case 0 :
						dx = (cutoff - point.x) * (cutoff - point.x);
						break;
					case 1 :
						dy = (cutoff - point.y) * (cutoff - point.y);
						break;
					case 2 :
						dz = (cutoff - point.z) * (cutoff - point.z);
						break;
				}
				d = dx + dy + dz;
				if ((d < distance) && (firstChild != null)) {
					// look only if child is closer than radius
					p2 = firstChild.findNearest(point, distance, dx, dy, dz);
				}
				if (p2 >= 0)
					p = p2;
			}
			// get dimension
			// get dimensional split
			// get furthest subnode
			// get distances
		}
		return p;
	}

	// quicksort the subarraypointsleft .. right] using 3-way partitioning
	private void sortKH(int left, int right, int dim) {

		// System.out.println("from KD Node (sorting from "+left+" to "+right+" by
		// "+dim+")" + "right<=left: " + (right <= left));
		if (right <= left)
			return;
		int lt = left, gt = right;
		Point3d temp = points[left];
		int i = left;
		int cmp = 0;
		while (i <= gt) {
			switch (dim) {
				case 0 :
					cmp = compareTo_X(points[i], temp);
					break;
				case 1 :
					cmp = compareTo_Y(points[i], temp);
					break;
				case 2 :
					cmp = compareTo_Z(points[i], temp);
			}

			if (cmp < 0)
				exch(points, lt++, i++);
			else if (cmp > 0)
				exch(points, i, gt--);
			else
				i++;
		}

		// points[left..lt-1] < temp = points[lt..gt] < points[gt+1..right].
		sortKH(left, lt - 1, dim);
		sortKH(gt + 1, right, dim);
		switch (dim) {
			case 0 :
				assert isSortedX(points, left, right);
				break;
			case 1 :
				assert isSortedY(points, left, right);
				break;
			case 2 :
				assert isSortedZ(points, left, right);
		}
	}

	/***********************************************************************
	 * Helper sorting functions
	 ***********************************************************************/
	// compareTo_...(Point3d v, Point3d w)
	// replaces:
	// ---> points[i].compareTo(temp)
	// ---> v.compareTo(w)

	private int compareTo_X(Point3d v, Point3d w) {
		int result;
		if (v.x < w.x) {
			result = -1;
			return result;
		} else if (v.x == w.x) {
			result = 0;
			return result;
		} else {
			result = 1;
			return result;
		}
	}

	private int compareTo_Y(Point3d v, Point3d w) {
		int result;
		if (v.y < w.y) {
			result = -1;
			return result;
		} else if (v.y == w.y) {
			result = 0;
			return result;
		} else {
			result = 1;
			return result;
		}
	}

	private int compareTo_Z(Point3d v, Point3d w) {
		int result;
		if (v.z < w.z) {
			result = -1;
			return result;
		} else if (v.z == w.z) {
			result = 0;
			return result;
		} else {
			result = 1;
			return result;
		}
	}

	// exchange points[i] and points[j]
	private void exch(Point3d[] points, int i, int j) {
		Point3d swap = points[i];
		points[i] = points[j];
		points[j] = swap;
	}

	/***********************************************************************
	 * Check if array is sorted - useful for debugging
	 ***********************************************************************/

	private boolean isSortedX(Point3d[] points, int left, int right) {
		for (int i = left + 1; i <= right; i++) {
			if (points[i].x < points[i - 1].x)
				return false;
		}
		return true;
	}

	private boolean isSortedY(Point3d[] points, int left, int right) {
		for (int i = left + 1; i <= right; i++) {
			if (points[i].y < points[i - 1].y)
				return false;
		}
		return true;
	}

}
