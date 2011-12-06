public class AACluster {
	
	public static boolean VVV = false;
	public static boolean VV = false;
        public static boolean V = false;
	private static Models[] m =  Models.values();
	
	/**
	 * @param args none
	 */
	public static void main(String[] args) {
		double[][] d1, d2, d3, d4;
		AAModel a, b;
		
		if(VV) {
			for(Models d : Models.values()) {
				AAModel tmp = new AAModel(d);
				System.out.println(tmp.toString());
			}
			System.exit(0);
		} else if (V) {
			System.out.println(new AAModel(Models.MTMAM).toString());
			System.exit(0);
		}

		d1 = new double[m.length][m.length];
		for(int i = 0; i < m.length; i++) {
			a = new AAModel(m[i]).scaleMax();
			for(int j = 0; j < m.length; j++) {
				b = new AAModel(m[j]).scaleMax();
				d1[i][j] = AAModel.dist(a, b);
			}
		}
		
		d2 = new double[m.length][m.length];
		for(int i = 0; i < m.length; i++) {
			a = new AAModel(m[i]).scaleOcc(20);
			for(int j = 0; j < m.length; j++) {
				b = new AAModel(m[j]).scaleOcc(20);
				d2[i][j] = AAModel.dist(a, b);
			}
		}
		
		d3 = new double[m.length][m.length];
		for(int i = 0; i < m.length; i++) {
			a = new AAModel(m[i]).scaleOne();
			for(int j = 0; j < m.length; j++) {
				b = new AAModel(m[j]).scaleOne();
				d3[i][j] = AAModel.dist(a, b);
			}
		}
	
//		d4 = new double[m.length][m.length];
//               for(int i = 0; i < m.length; i++) {
//                        a = new AAModel(m[i]).scaleOneFLess(1);
//                        for(int j = 0; j < m.length; j++) {
//                                b = new AAModel(m[j]).scaleOneFLess(1);
//                               d4[i][j] = AAModel.dist(a, b);
//                        }
//                }
	
		System.out.println("Scaled by maximum: ");
		printDistMat(d1);
		System.out.println("Scaled by number of occurence of AAs: ");
		printDistMat(d2);
		System.out.println("Scaled to one subst. per time step: ");
		printDistMat(d3);
//		System.out.println("Scaled to one disregarding AA Frequencies: ");
//		printDistMat(d4);
		
//		System.out.println("DAYHOFF - DCMUT: " + AAModel.dist(new AAModel("DAYHOFF").scaleMax(), new AAModel("DCMUT").scaleMax()));
//		System.out.println("JTT - JTTDCMUT: " + AAModel.dist(new AAModel("JTT").scaleMax(), new AAModel("JTTDCMUT").scaleMax()));
//		System.out.println("MTART - MTZOA: " + AAModel.dist(new AAModel("MTART").scaleMax(), new AAModel("MTZOA").scaleMax()));
	}

	private static void printDistMat(double[][] d) {
		String s = "     ";
		
		for(int i = 0; i < m.length; i++) {
			int l = m[i].name().length() > 4 ? 4 : m[i].name().length();
			s += String.format("%4s ", m[i].name().substring(0, l));
		}
		
		s+= "\n";
		
		for(int i = 0; i < d.length; i++) {
			int l = m[i].name().length() > 4 ? 4 : m[i].name().length();
			s += String.format("%4s ", m[i].name().subSequence(0, l));
			for(int j = 0; j < d[i].length; j++) {
				s += String.format("%4.0f ", d[i][j]);
			}
			s += "\n";
		}
		System.out.println(s);
	}
}
