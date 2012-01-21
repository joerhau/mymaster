import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

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
		
		for(Models d : Models.values()) {
			FileWriter fstream;
			AAModel tmp = new AAModel(d);
			tmp.scaleMax();
			tmp.revertMax();
//			System.out.println(tmp.matToString());
			try {
				fstream = new FileWriter(d.name() + ".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				out.write(tmp.matToString());
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
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
	
		d4 = new double[m.length][m.length];
		for(int i = 0; i < m.length; i++) {
			a = new AAModel(m[i]).scaleOneFLess();
			for(int j = 0; j < m.length; j++) {
				b = new AAModel(m[j]).scaleOneFLess();
				d4[i][j] = AAModel.dist(a, b);
			}
		}
		
//		AAModel[] v = new AAModel[3];
//		v[0] = new AAModel(Models.FLU).scaleMax();
//		v[1]  = new AAModel(Models.HIVB).scaleMax();
//		v[2] = new AAModel(Models.HIVW).scaleMax();
//
//		
//		AAModel[] d = new AAModel[4];
//		d[0] = new AAModel(Models.DAYHOFF).scaleMax();
//		d[1]  = new AAModel(Models.DCMUT).scaleMax();
//		d[2] = new AAModel(Models.JTT).scaleMax();
//		d[3] = new AAModel(Models.JTTDCMUT).scaleMax();
//		
//		AAModel viral = AAModel.getRepresentative(v, "Vir").scaleMax();
//		AAModel dayhoffs = AAModel.getRepresentative(d, "Day").scaleMax();
//		
//		System.out.println(viral.toString());
//		System.out.println(dayhoffs.toString());
//		AAModel tmp = new AAModel(Models.WAG).scaleMax();
		
//		for(int i = 0; i < v.length ; i++) {
//			System.out.println(v[i].name + " - " + viral.name + ": \t" + AAModel.dist(v[i], viral));
//			System.out.println(viral.name + " - " + tmp.name + ": \t" + AAModel.dist(viral, tmp));
//		}
		
//		for(int i = 0; i < d.length ; i++) {
//			System.out.println(d[i].name + " - " + dayhoffs.name + ": \t" + AAModel.dist(d[i], dayhoffs));
//			System.out.println(dayhoffs.name + " - " + tmp.name + ": \t" + AAModel.dist(dayhoffs, tmp));
//		}
		
		BufferedWriter max, aac, one, oneF;
		

		
//		System.out.println("Scaled by maximum: ");
//		d1 = scaleByMax(d1, 100);
//		System.out.println("Scaled by number of occurence of AAs: ");
//		d2 = scaleByMax(d2, 100);
//		distMatToString(d2);
//		System.out.println("Scaled to one subst. per time step: ");
//		d3 = scaleByMax(d3, 100);
//		distMatToString(d3);
//		System.out.println("F less: ");
//		d4 = scaleByMax(d4, 100);
//		distMatToString(d4);
		
//		try{
//			max = new BufferedWriter(new FileWriter("max.txt"));
//			aac = new BufferedWriter(new FileWriter("aac.txt"));
//			one = new BufferedWriter(new FileWriter("one.txt"));
//			oneF = new BufferedWriter(new FileWriter("oneF.txt"));
//			max.write(distMatToString(d1));
//			aac.write(distMatToString(d2));
//			one.write(distMatToString(d3));
//			oneF.write(distMatToString(d4));
//			max.close();
//			aac.close();
//			one.close();
//			oneF.close();
//		}catch (Exception e){//Catch exception if any
//			System.err.println("Error: " + e.getMessage());
//		}
		
//		System.out.println(new AAModel(Models.FLU).toString());
//		System.out.println(new AAModel(Models.HIVB).toString());
//		System.out.println(new AAModel(Models.LG).toString());
		
//		System.out.println("DAYHOFF - DCMUT: " + AAModel.dist(new AAModel("DAYHOFF").scaleMax(), new AAModel("DCMUT").scaleMax()));
//		System.out.println("JTT - JTTDCMUT: " + AAModel.dist(new AAModel("JTT").scaleMax(), new AAModel("JTTDCMUT").scaleMax()));
//		System.out.println("MTART - MTZOA: " + AAModel.dist(new AAModel("MTART").scaleMax(), new AAModel("MTZOA").scaleMax()));
	}
	
	private static String distMatToString(double[][] d) {
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
		return s;
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
