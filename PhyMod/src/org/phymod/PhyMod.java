package org.phymod;

import java.util.Random;

import org.phymod.io.PartitionLoader;

public class PhyMod {
	private static PhylipAlignment phy;
	private static PartitionLoader part;
	public static boolean del = false;

	public static void main(String[] args) {
		int count = 0;
		String outfile = "";
		
		if(args.length < 2) { 
			System.out.println("Usage: PhyMod [Options] PHYLIPFILE PARTITIONFILE");
			System.exit(1);
		}
		
		part = new PartitionLoader(args[1]);
		phy = new PhylipAlignment(args[0], part);
		
		for(int i = 2; i < args.length; i++) {
			if(args[i].substring(0,1).equals("--"));
			else if(args[i].substring(0,1).equals("-")) {
				if(args[i].substring(1,2).equals("m")) count = Integer.parseInt(args[++i]);
				else if(args[i].substring(1,2).equals("n")) outfile = args[++i];
				else if(args[i].substring(1,2).equals("d")) del = true;
			}
		}
		
		PhylipAlignment rand = randomPartitions(count);
		
		if(!outfile.equals("")) {
			rand.toFile(outfile);
		}
		else System.out.println(rand.phylipToString());
	}
	
	public static PhylipAlignment randomPartitions(int count) {
	    Species[] s = new Species[phy.spec.length];
	    int[] r = createRand(count);
	    String str = "";
	    for(int i = 0; i < r.length; i++)
	    	str += r[i] + " ";
	    
	    System.out.println("Randomly chosen partitions " + str);
	    
	    for (int i = 0; i < phy.spec.length; i++) {
			s[i] = new Species(phy.spec[i].name, count);
		}
			
		for (int i = 0; i < phy.spec.length; i++) {
			for(int j = 0; j < r.length; j++) {
				s[i].replace(phy.spec[i], r[j], j);
			}
		}
		
		return new PhylipAlignment(s);
	}
	
	private static int[] createRand(int count) {
		Random rand = new Random();
		boolean redo;
		int[] r = new int[count];
		
		do {
			redo = false;
			for(int i = 0; i < count; i++) {
				boolean rep;
		    	int tmp;
		    	
		    	do {
		    		rep = false;
		    		tmp = rand.nextInt(phy.nrPartitions);
		    		for(int j = 0; j < i; j++) {
		    			if(r[j] == tmp) {
		    				rep = true;
		    				break;
		    			}
		    		}
		    	} while(rep);
		    	
		    	r[i] = tmp;
		    }
		
			for(int j = 0; j < phy.spec.length; j++) {
				String str = "";
				
				for(int i = 0; i < r.length; i++)
					str += phy.spec[j].partitions[r[i]].data;
				
				if(str.replaceAll("-", "").equals("")) {
					redo = true;
					break;
				}
			}
		} while (redo);
		
		return r;
	}
}
