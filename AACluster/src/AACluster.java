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
		
		for(Models d : Models.values())
			System.out.println(new AAModel(d).toString());
		
		DecimalFormat dd = new DecimalFormat("#0000");
		String s = "     ";
		
		for(int i = 0; i < m.length; i++) {
			int l = m[i].name().length() > 4 ? 4 : m[i].name().length();
			s += String.format("%-4s ", m[i].name().substring(0, l));
		}
		
		s+= "\n";
		
		for(int i = 0; i < dists.length; i++) {
			int l = m[i].name().length() > 4 ? 4 : m[i].name().length();
			s += String.format("%-4s ", m[i].name().subSequence(0, l));
			for(int j = 0; j < dists[i].length; j++) {
//				s += dists[i][j] + " ";
				s+= dd.format(dists[i][j]) + " ";
			}
			s += "\n";
		}
		System.out.println(s);
		
		System.out.println("DAYHOFF - DCMUT: " + AAModel.dist(new AAModel("DAYHOFF"), new AAModel("DCMUT")));
		System.out.println("JTT - JTTDCMUT: " + AAModel.dist(new AAModel("JTT"), new AAModel("JTTDCMUT")));
		System.out.println("MTART - MTZOA: " + AAModel.dist(new AAModel("MTART"), new AAModel("MTZOA")));
	}

}
