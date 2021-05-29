import java.util.Random;
import java.io.BufferedReader;
import java.io.FileReader;

/**
 * Do very rudimentary timing of the getNearest method, as well as very small sampling to measure
 * accuracy. You may wish to increase the amount of both.
 *
 * @author brinkmwj
 * @version 2009-10-15
 */

public class Controller {
	/**
	 * This is a helper function that reads in the Starbucks.csv file. If you modified your Starbucks.csv,
	 *  or if Starbucks.csv is not in the project's folder, then this might blow up on you.
	 */
	public static Starbucks.StarbucksLocation[] readEntryList(){
		String fname = "starbucks.csv";
		Starbucks.StarbucksLocation[] ret = null;

		try {
			BufferedReader br = new BufferedReader(new FileReader(fname));
			String line;
			line = br.readLine(); //First line is the count of number of Starbucks locations
			String[] components = line.split(",");

			int count = Integer.parseInt(components[0]);
			ret = new Starbucks.StarbucksLocation[count];

			line = br.readLine(); //Column labels, to be discarded
			for(int i=0; i < count; i++){
				line = br.readLine();
				components = line.split(",");
				ret[i] = new Starbucks.StarbucksLocation(components[0], components[3], 
						Double.parseDouble(components[2]), Double.parseDouble(components[1]));
			}
		} catch (Exception e) {
			e.printStackTrace(System.out);
		}
		return ret;
	}

	public static Starbucks.StarbucksLocation bruteGetNearest(double lng, double lat, Starbucks.StarbucksLocation[] slocs){
		double dist = Starbucks.distance(slocs[0].lng, slocs[0].lat, lng, lat);
		int index = 0;
		for(int i=1; i < slocs.length; i++){
			double cdist = Starbucks.distance(slocs[i].lng, slocs[i].lat, lng, lat);
			if(cdist < dist){
				index = i;
				dist = cdist;
			}
		}
		return new Starbucks.StarbucksLocation(slocs[index]);
	}

	public static void main(String[] args){
		Starbucks sS = new StudentStarbucks();

		Starbucks.StarbucksLocation[] slocs = readEntryList();
		Starbucks.StarbucksLocation[] slocsCopy = new Starbucks.StarbucksLocation[slocs.length];

		Random rng = new Random(5);
		//randomly scramble the entryList
		for(int i=0; i<slocs.length; i++){
			int j = i + (int)(rng.nextDouble()*(slocs.length-i));
			Starbucks.StarbucksLocation t = slocs[i];
			slocs[i] = slocs[j];
			slocs[j] = t;
		}
		//Copy the list, so I have a clean copy for accuracy testing
		for(int i=0; i < slocs.length; i++){
			slocsCopy[i] = new Starbucks.StarbucksLocation(slocs[i]);
		}

		/*
		 * Unless your build() method is very slow this time will not be very accurate.
		 */
		long startb = System.nanoTime();
		sS.build(slocs);
		long endb = System.nanoTime();
		System.out.println("Building the data structure took: " + (endb-startb)/1000000.0 + " milliseconds (time not very accurate)");

		/*
		 * Quick accuracy test for some known values - Only use if your code is supposed to be perfectly accurate!
                 */
		Starbucks.StarbucksLocation s = sS.getNearest(-76.16304, 39.5106);
		if(!s.address.equals("Maryland House Travel Plaza")){
			System.out.println("Failed Maryland House Travel Plaza test. Got: " + s.address);
			System.exit(0);
		}
    		s = sS.getNearest(-104.8287, 39.70294);
		if(!s.address.equals("I-225 - Alameda-Aurora")){
			System.out.println("Failed I-225 - Alameda-Aurora test. Got: " + s.address);
			System.exit(0);
		}
		
		/* 
		 * TEST FOR SPEED
		 * Do many searches, and compute the average time per search. -- This should take between
		 * 5 and 50 seconds.
		 *
 		 * On my machine with my solutions I got the following time/accuracy:
		 * Brute force: 0.162 ms per search, 1.08 ms build, 0.0 error rate
		 * My solution: 0.026 ms per search, 39.3 ms build, 0.0 error rate
		 */
		long start = 0;
		long end = 0;
		int numTrials = 1000;
		do{
			start = System.nanoTime();
			for(int i=0; i<numTrials; i++){
				double x = -125.0 + 73.0*rng.nextDouble();
				double y = 24.0 + 25.0*rng.nextDouble();
				Starbucks.StarbucksLocation tmp = sS.getNearest(x, y);
//				Starbucks.StarbucksLocation opt = bruteGetNearest(x, y, slocsCopy);

			}
			end = System.nanoTime();
			numTrials *= 10;
		} while ((end - start)/1000000.0 < 5000 && numTrials < 1000000);

		System.out.println("Time: " + (((end - start)/1000000.0)/numTrials) + " ms per search, " + numTrials + " trials");

		/* 
		 * TEST FOR ACCURACY
		 * Use brute-force for comparison
		 */
		double studentTotal = 0.0;
		double optTotal = 0.0;
		numTrials = 1000;
		for(int i=0; i < numTrials; i++){
			double x = -125.0 + 73.0*rng.nextDouble();
			double y = 24.0 + 25.0*rng.nextDouble();
			Starbucks.StarbucksLocation student = sS.getNearest(x, y);
			Starbucks.StarbucksLocation opt = bruteGetNearest(x, y, slocsCopy);
			studentTotal += Starbucks.distance(student.lng, student.lat, x, y);
			optTotal += Starbucks.distance(opt.lng, opt.lat, x, y);
		}

		double error = studentTotal/optTotal;
		System.out.println("Error percentage is: " + 100.0*(error-1.0));
	}
}
