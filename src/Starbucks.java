/*****
 * Author   : brinkmwj
 * Date     : 2016-11-07
 * Sources  : All code is original
 * Purpose  : This class is designed to be an "Interface." Your starbucks data structure should be
 *            a sub-class of this class, but you should not implement this class. Instead, you should
 *            create a StudentStarbucks class in the file StudentStarbucks.java
 * Note     : DO NOT CHANGE THIS FILE!!! 
 */

public abstract class Starbucks {

	public static class StarbucksLocation {
		public String city;
		public String address;
		public double lng; //longitude
		public double lat; //latitude

		public StarbucksLocation(){
			city = "";
			address = "";
			lng = lat = 0.0;
		}

		public StarbucksLocation(String icity, String iaddress, double ilng, double ilat){
			city = icity;
			address = iaddress;
			lng = ilng;
			lat = ilat;
		}

		public StarbucksLocation(StarbucksLocation orig){
			city = orig.city;
			address = orig.address;
			lng = orig.lng;
			lat = orig.lat;
		}
	}

	/**
	 * Given the longitude and latitude of two points, return their distance in km
	 *
	 * @param long1 longitude of first point
	 * @param lat1 latitude of first point
	 * @param long2 longitude of second point
	 * @param lat2 latitude of second point
	 * @return distance in km
	 */
	public static double distance(double long1, double lat1, double long2, double lat2){
		double R = 6371; // radius of earth, in km
		double dLat = (lat2-lat1)*Math.PI/180;
		double dLon = (long2-long1)*Math.PI/180;
		lat1 = lat1*Math.PI/180;
		lat2 = lat2*Math.PI/180;

		double a = Math.sin(dLat/2) * Math.sin(dLat/2) +
			Math.sin(dLon/2) * Math.sin(dLon/2) * Math.cos(lat1) * Math.cos(lat2); 
		double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a)); 
		double d = R * c;

		return d;
	}

	/**
	 * Add all entries in the array to your data structure
	 *
	 * Note: If you detect that two items have the same coordinates, then do not add the new item
	 *       that has the same coordinates as another item. This is guaranteed to happen, by the way,
	 *       because some Starbucks locations are listed in the database twice. We will define two locations
	 *       to be the "same location" if both |x1 - x2| <= 0.00001 and |y1 - y2| < 0.00001
	 *
	 * @param allLocations An array containing all of the Starbucks locations
	 */
	public abstract void build(StarbucksLocation[] allLocations);

	/*
	 * Return a reference to a deep copy of the entry that is closest to the given coordinates. Your
	 *  answer may be approximate, but then you will lose points on the "Accuracy" quality measure
	 */
	public abstract StarbucksLocation getNearest(double lng, double lat);
}
