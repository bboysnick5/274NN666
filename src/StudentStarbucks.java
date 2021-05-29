
public class StudentStarbucks extends Starbucks {

	
	private KDTree.Point transLatLngToXYZPt(double lng, double lat) {
		return new KDTree.Point(new double[]{
				Math.cos(Math.toRadians(lat)) * Math.cos(Math.toRadians(lng)),
				Math.cos(Math.toRadians(lat)) * Math.sin(Math.toRadians(lng)),
				Math.sin(Math.toRadians(lat))});
	}
	
	@Override
	public void build(StarbucksLocation[] allLocations) {
		for (StarbucksLocation loc : allLocations) {
			tree.insert(transLatLngToXYZPt(loc.lng, loc.lat), loc);
		}
		System.out.println("size " + allLocations.length);
	}

	@Override
	public StarbucksLocation getNearest(double lng, double lat) {
		return tree.nearestValue(transLatLngToXYZPt(lng, lat));
	}

	
	public StudentStarbucks() {
		super();
		tree = new KDTree<StarbucksLocation>(3);
	}

	private KDTree<StarbucksLocation> tree;	
	
}
