import java.util.Arrays;

public class KDTree<V> {
	
	public KDTree(int dims) {
		this.DIMS = dims;
	}
	
	
	public static class Point {
		public double[] coords;
		
		public Point(double[] coords) {
			this.coords = new double[coords.length];
			this.coords = Arrays.copyOf(coords, coords.length);
		}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj) 
				return true;
			if (!(obj instanceof Point)) 
				return false;
			Point otherPoint = (Point)obj;
			return Arrays.equals(this.coords, otherPoint.coords);
		}
		
		public static double EucDist(Point p1, Point p2) {
			double result = 0.0;
		    for (int i = 0; i < p1.coords.length; ++i)
		        result += (p1.coords[i] - p2.coords[i]) * (p1.coords[i] - p2.coords[i]);
		    return Math.sqrt(result);
		}
	}
	
	private class Node<V> {
		public Point p;
		public V value;
		public Node<V> left;
		public Node<V> right;
		
		public Node(Point p, V v) {
			this.p = p;
			this.value = v;
			left = right = null;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj) 
				return true;
			if (!(obj instanceof Node)) 
				return false;
			Node<V> otherNode = (Node<V>)obj;
			return this.p.equals(otherNode.p);
		}
		
	}
	
	public void insert(Point p, V v) {
		if (root == null) {
			root = new Node<V>(p, v);
			return;
		}
		Node<V> cur = root;
		Node<V> parent = null;
		int level = 0;
		while (true) {
			if (cur.p.equals(p)) {
				cur.value = v;
				return;
			}
			parent = cur;
			if (p.coords[level%DIMS] < cur.p.coords[level%DIMS]) {
				cur = cur.left;
				if (cur == null) {
					parent.left = new Node<V>(p, v);
					return;
				}
			} else {
				cur = cur.right;
				if (cur == null) {
					parent.right = new Node<V>(p, v);
					return;
				}
			}
			level++;
		}
	}
	
	public boolean contains(Point p) {
		return findValue(p) != null;
	}
	
	public V findValue(Point p) {
		Node<V> cur = root;
		int level = 0;
		while(cur != null) {
			if (cur.p.equals(p)) {
				return cur.value;
			}
			cur = p.coords[level%DIMS] < cur.p.coords[level%DIMS] ? cur.left : cur.right; 
		}
		return null;
	}
	
	private void nearestValueHelper(Node<V> cur, Point p, int level) {
		if (cur == null)
			return;
		double curDist = Point.EucDist(cur.p, p);
		if (curDist < bestDist) {
			bestDist = curDist;
			guess = cur;
		}
		Node<V> next = p.coords[level%DIMS] < cur.p.coords[level%DIMS] ? cur.left : cur.right;
		Node<V> other = next == cur.left ? cur.right : cur.left;
		nearestValueHelper(next, p, level+1);
		if (Math.abs(p.coords[level%DIMS] - cur.p.coords[level%DIMS]) < bestDist)
			nearestValueHelper(other, p, level+1);
		return;
	}
	
	public V nearestValue(Point p) {		
		guess = null;
		bestDist = Double.MAX_VALUE;
		nearestValueHelper(root, p, 0);
		return guess.value;
	}
	
	private Node<V> guess;
	private double bestDist;
	
	private Node<V> root;
	private final Integer DIMS;
	
}
