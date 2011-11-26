import java.text.DecimalFormat;


public class AACluster {

	/**
	 * @param args none
	 */
	public static void main(String[] args) {
		double[][] dists;
		
		AAModel a, b;
		Models[] m = Models.values();
		dists = new double[m.length][m.length];
		for(int i = 0; i < m.length; i++) {
			a = new AAModel(m[i].name());
			for(int j = 0; j < m.length; j++) {
				b = new AAModel(m[j].name());
				dists[i][j] = AAModel.dist(a, b);
			}
		}
		
		DecimalFormat dd = new DecimalFormat("#00");
		String s = "";
		for(int i = 0; i < dists.length; i++) {
			for(int j = 0; j < dists[i].length; j++) 
					s+= dd.format(dists[i][j]) + " ";
			s += "\n";
		}
		System.out.println(s);
		
		System.out.println("DAYHOFF - DCMUT: " + AAModel.dist(new AAModel("DAYHOFF"), new AAModel("DCMUT")));
		System.out.println("JTT - JTTDCMUT: " + AAModel.dist(new AAModel("JTT"), new AAModel("JTTDCMUT")));
	}

}
