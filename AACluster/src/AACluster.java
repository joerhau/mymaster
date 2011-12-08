public class AACluster {
	
	public static boolean VVV = false;
	public static boolean VV = false;
        public static boolean V = false;
	private static Models[] m =  Models.values();
	
	/**
	 * @param args none
	 */
	public static void main(String[] args) {
		double[][] d1, d2, d3;
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
	
		System.out.println("Scaled by maximum: ");
		d1 = scaleByMax(d1, 100);
		printDistMat(d1);
		System.out.println("Scaled by number of occurence of AAs: ");
		d2 = scaleByMax(d2, 100);
		printDistMat(d2);
		System.out.println("Scaled to one subst. per time step: ");
		d3 = scaleByMax(d3, 100);
		printDistMat(d3);
		
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
	
	/**
	 * scales all entries in mat so that there are only values between 0 and 1
	 * 
	 * @param mat
	 * @return
	 */
	public static double[][] scaleByMax(double[][] mat, double scaler) {
		double max = 0;
		
		for (int i = 0; i < mat.length; i++)
			for (int j = 0; j < mat.length; j++)
				if (mat[i][j] > max)
					max = mat[i][j];
		
		double s =  scaler / max;
		
		for (int i = 0; i < mat.length; i++)
			for (int j = 0; j < mat.length; j++)
				mat[i][j] *= s;
		return mat;
	}
	
	public static double dist(double[][] a, double[][] b) {
		double ssum = 0;
		for(int i = 0; i < a.length; i++) {
			for(int j = i + 1; j < a[i].length; j++) {
				ssum += Math.pow(a[i][j] - b[i][j], 2);
			}
		}
		return Math.sqrt(ssum);
	}
	

}
